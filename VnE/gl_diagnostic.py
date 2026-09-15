"""Opt-in GL diagnostics, usable unchanged with both ASID GL implementations.

Run from any directory: python VnE/gl_diagnostic.py
Close the windows when finished. No backend or renderer settings are changed.
"""
import argparse
import ctypes as C
import ctypes.util
import functools
import json
import os
from pathlib import Path
import sys
import traceback

import PySide6
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtOpenGL import QOpenGLVersionFunctionsFactory, QOpenGLVersionProfile

ROOT = Path(__file__).resolve().parent
SEEN = set()
ACTIVE = []
COUNTS = {"compile": 0, "link": 0, "failed": 0}
WRAPPER = None
WIDGETS = []
FAILED_WIDGETS = set()


def emit(event, **data):
    print("GLDIAG " + json.dumps(dict(event=event, **data), default=str), flush=True)


def once(key, event, **data):
    if key not in SEEN:
        SEEN.add(key)
        emit(event, **data)


def identity(obj):
    return None if obj is None else hex(id(obj))


def fmt(value):
    fields = ("majorVersion", "minorVersion", "profile", "renderableType",
              "swapBehavior", "samples", "depthBufferSize", "stencilBufferSize",
              "alphaBufferSize", "options")
    return {name: getattr(value, name)() for name in fields}


def proc(name, result, *args):
    # Resolve against the *current* Qt context, without PyOpenGL or ASID caches.
    ctx = QtGui.QOpenGLContext.currentContext()
    if ctx is None:
        raise RuntimeError("No current Qt context for " + name)
    address = ctx.getProcAddress(name.encode())
    if not address:
        raise RuntimeError("Qt cannot resolve " + name)
    factory = C.WINFUNCTYPE if os.name == "nt" else C.CFUNCTYPE
    return factory(result, *args)(int(address))


def integer(token, count=1):
    values = (C.c_int * count)()
    proc("glGetIntegerv", None, C.c_uint, C.POINTER(C.c_int))(token, values)
    return list(values) if count != 1 else values[0]


def errors():
    get = proc("glGetError", C.c_uint)
    values = []
    for _ in range(16):
        value = get()
        if not value:
            break
        values.append(hex(value))
    return values


def native_contexts():
    result = {}
    if not sys.platform.startswith("linux"):
        return {"native": "EGL/GLX probes only on Linux"}
    for library, name in (("EGL", "eglGetCurrentContext"), ("GLX", "glXGetCurrentContext")):
        try:
            path = C.util.find_library(library)
            if not path:
                result[library] = "library unavailable"
                continue
            lib = C.CDLL(path)
            function = getattr(lib, name)
            function.restype, function.argtypes = C.c_void_p, []
            result[library] = hex(function() or 0)
        except Exception as exc:
            result[library] = str(exc)
    return result


def backend(widget):
    ctx = widget.context()
    if ctx is None:
        return
    key = ("backend", id(widget), id(ctx))
    if key in SEEN:
        return
    SEEN.add(key)
    data = dict(widget=type(widget).__name__, context=identity(ctx), valid=ctx.isValid(),
                isOpenGLES=ctx.isOpenGLES(), actual=fmt(ctx.format()),
                widgetFormat=fmt(widget.format()), native=native_contexts())
    get = proc("glGetString", C.c_char_p, C.c_uint)
    data["GL"] = {name: (get(token) or b"").decode(errors="replace") for name, token in
                  (("vendor", 0x1F00), ("renderer", 0x1F01), ("version", 0x1F02),
                   ("GLSL", 0x8B8C))}
    try:
        functions = QOpenGLVersionFunctionsFactory.get(QOpenGLVersionProfile(ctx.format()), ctx)
        data["QtFunctionTable"] = type(functions).__name__
    except Exception as exc:
        data["QtFunctionTable"] = {"error": str(exc)}
    # Inspect PyOpenGL only if the application already imported it. Importing it
    # in the Qt-only version would change the very loader path being measured.
    platform = sys.modules.get("OpenGL.platform")
    if platform:
        p = platform.PLATFORM
        try:
            data["PyOpenGL"] = dict(version=sys.modules["OpenGL"].__version__,
                                    platform=type(p).__module__ + "." + type(p).__name__,
                                    current=hex(p.GetCurrentContext() or 0),
                                    GL_library=getattr(p.GL, "_name", None),
                                    contextLookup=str(p.GetCurrentContext),
                                    resolver=str(p.getExtensionProcedure))
        except Exception as exc:
            data["PyOpenGL"] = {"error": str(exc)}
    if WRAPPER:
        selected = WRAPPER.openglFunctions
        data["wrapper"] = dict(selected=identity(selected),
                               owner=identity(selected.context.context()) if selected is not None else None,
                               contexts=len(WRAPPER.contexts))
    emit("backend", **data)


def check(stage, widget=None):
    ctx = QtGui.QOpenGLContext.currentContext()
    widget = widget or (ACTIVE[-1] if ACTIVE else None)
    data = dict(stage=stage, widget=type(widget).__name__ if widget else None,
                current=identity(ctx), expected=identity(widget.context()) if widget else None,
                currentMatches=(ctx is widget.context()) if widget else None,
                threadMatches=(ctx.thread() is QtCore.QThread.currentThread()) if ctx else None)
    if ctx:
        data.update(drawFBO=integer(0x8CA6), readFBO=integer(0x8CAA),
                    widgetFBO=widget.defaultFramebufferObject() if widget else None,
                    viewport=integer(0x0BA2, 4), samples=integer(0x80A9),
                    drawBuffer=hex(integer(0x0C01)), readBuffer=hex(integer(0x0C02)),
                    fboStatus=hex(proc("glCheckFramebufferStatus", C.c_uint, C.c_uint)(0x8D40)),
                    errors=errors(), shaders=dict(COUNTS))
        if WRAPPER:
            selected = WRAPPER.openglFunctions
            data["wrapperOwner"] = identity(selected.context.context()) if selected is not None else None
            data["wrapperMatches"] = selected is not None and selected.context.context() is ctx
    emit("stage", **data)


def shader_status(name, object_id):
    shader = name == "glCompileShader"
    kind = "Shader" if shader else "Program"
    value = C.c_int()
    query = proc("glGet" + kind + "iv", None, C.c_uint, C.c_uint, C.POINTER(C.c_int))
    query(object_id, 0x8B81 if shader else 0x8B82, C.byref(value))
    ok = bool(value.value)
    query(object_id, 0x8B84, C.byref(value))
    log = C.create_string_buffer(max(1, value.value))
    proc("glGet" + kind + "InfoLog", None, C.c_uint, C.c_int,
         C.POINTER(C.c_int), C.c_void_p)(object_id, len(log), None, log)
    COUNTS["compile" if shader else "link"] += 1
    COUNTS["failed"] += int(not ok)
    if not ok or log.value:
        emit("shader", kind=kind, object=int(object_id), ok=ok,
             log=log.value.decode(errors="replace"))


def trace_gl(name, function, owner=None):
    @functools.wraps(function)
    def call(*args, **kwargs):
        ctx = QtGui.QOpenGLContext.currentContext()
        expected = ACTIVE[-1].context() if ACTIVE else None
        mismatch = ctx is None or (expected is not None and ctx is not expected)
        if owner is not None and ctx is not owner.context.context():
            mismatch = True
        if mismatch:
            once(("mismatch", name, identity(ctx), identity(expected)), "GL-context-mismatch",
                 function=name, current=identity(ctx), expected=identity(expected),
                 wrapperOwner=identity(owner.context.context()) if owner else None)
        result = function(*args, **kwargs)
        if name in ("glCompileShader", "glLinkProgram") and ctx is not None:
            shader_status(name, args[0])
        return result
    return call


def instrument(main):
    global WRAPPER
    WRAPPER = sys.modules.get("Source.opengl")
    if WRAPPER:
        original = WRAPPER.OpenGLFunctions.__getattr__

        def lookup(self, name):
            result = original(self, name)
            return trace_gl(name, result, self) if name.startswith("gl") and callable(result) else result

        WRAPPER.OpenGLFunctions.__getattr__ = lookup
    else:
        for module_name, module in list(sys.modules.items()):
            if module_name == "Source" or module_name.startswith("Source."):
                for name, value in list(vars(module).items()):
                    if name.startswith("gl") and callable(value):
                        setattr(module, name, trace_gl(name, value))

    def wrap_method(cls, name, label):
        original = getattr(cls, name)

        @functools.wraps(original)
        def method(self, *args, **kwargs):
            widget = self if isinstance(self, QOpenGLWidget) else None
            if widget and id(widget) in FAILED_WIDGETS:
                return
            key = (id(self), label, identity(widget.context()) if widget else None)
            first = key not in SEEN
            SEEN.add(key)
            if widget:
                ACTIVE.append(widget)
                if widget not in WIDGETS:
                    WIDGETS.append(widget)
            try:
                if first:
                    if widget and widget.context() and QtGui.QOpenGLContext.currentContext() is widget.context():
                        backend(widget)
                    check(label + ":begin", widget)
                result = original(self, *args, **kwargs)
                if first:
                    check(label + ":end", widget)
                return result
            except BaseException as exc:
                once(("exception", label, str(exc)), "exception", stage=label, error=repr(exc))
                if widget and not isinstance(exc, KeyboardInterrupt):
                    # Keep the independent clear test alive even if ASID aborts
                    # initialization. Never render a partially initialized scene.
                    FAILED_WIDGETS.add(id(widget))
                    return
                raise
            finally:
                if widget:
                    ACTIVE.pop()
        setattr(cls, name, method)

    classes = [main.OpenGlWidget]
    drawer = sys.modules.get("Source.Extensions.ChemPack.DrawerWidget")
    if drawer:
        classes.append(drawer.DrawerGL)
    for cls in classes:
        for name in ("initializeGL", "resizeGL", "paintGL"):
            wrap_method(cls, name, cls.__name__ + "." + name)
    from Source import ShaderDataObjects, ShaderPipelines, UniformBuffers
    for cls, name in ((ShaderPipelines.ShaderProgramsCreator, "createProgram"),
                      (ShaderDataObjects.ShaderData, "__init__"),
                      (UniformBuffers.SceneUniformBuffer, "__init__")):
        # Shader creation is counted per shader above; summarize the stage once.
        wrap_method(cls, name, cls.__name__)
    wrap_method(ShaderPipelines.aShaderPipeline, "draw", "pipeline.draw")


class ClearProbe(QOpenGLWidget):
    def __init__(self, surface, loader=False):
        super().__init__()
        self.setFormat(surface)
        self.loader = loader
        self.functions = None
        self.failed = False
        self.frames = 0
        self.setWindowTitle("ASID magenta: " + ("application loader" if loader else "Qt getProcAddress"))
        self.resize(320, 220)
        self.frameSwapped.connect(self.swapped)

    def swapped(self):
        self.frames += 1

    def initializeGL(self):
        backend(self)
        if self.loader and WRAPPER:
            try:
                self.functions = WRAPPER.addContext(self)
            except Exception as exc:
                self.failed = True
                emit("clear-loader-failed", error=repr(exc))

    def makeCurrent(self, use_wrapper=True):
        if self.functions is not None and use_wrapper:
            return self.functions.makeCurrent()
        return super().makeCurrent()

    def paintGL(self):
        if self.failed:
            return
        # Exactly a clear, with no ASID initialization/shaders/buffers.
        if self.loader and WRAPPER:
            self.functions.glClearColor(1., 0., 1., 1.)
            self.functions.glClear(0x4000 | 0x0100)
        elif self.loader:
            from OpenGL.GL import glClearColor, glClear
            glClearColor(1., 0., 1., 1.)
            glClear(0x4000 | 0x0100)
        else:
            proc("glClearColor", None, C.c_float, C.c_float, C.c_float, C.c_float)(1., 0., 1., 1.)
            proc("glClear", None, C.c_uint)(0x4000 | 0x0100)
        key = ("clear", id(self))
        if key not in SEEN:
            SEEN.add(key)
            check("clear:" + ("application-loader" if self.loader else "Qt"), self)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clear-only", action="store_true", help="Qt-only clear, no ASID imports")
    parser.add_argument("--seconds", type=int, default=0, help="auto-close after N seconds (default: manual)")
    args = parser.parse_args()
    os.chdir(ROOT)  # VnE uses relative font and extension paths.
    sys.path.insert(0, str(ROOT))
    emit("environment", python=sys.version, PySide6=PySide6.__version__, Qt=QtCore.qVersion(),
         environment={k: os.environ.get(k) for k in
                      ("XDG_SESSION_TYPE", "WAYLAND_DISPLAY", "DISPLAY", "QT_QPA_PLATFORM",
                       "PYOPENGL_PLATFORM", "QT_XCB_GL_INTEGRATION", "QT_OPENGL", "LIBGL_ALWAYS_SOFTWARE")})
    app_main = None
    if not args.clear_only:
        try:
            from Source import MainWindow as app_main
        except Exception:
            emit("application-import-failed", traceback=traceback.format_exc())
    emit("before-QApplication", defaultFormat=fmt(QtGui.QSurfaceFormat.defaultFormat()))
    app = QtWidgets.QApplication([sys.argv[0]])
    emit("QApplication", platform=QtGui.QGuiApplication.platformName())
    window = None
    if app_main:
        window = app_main.MainWindow()
        instrument(app_main)
        surface = window.opengl_widget.format()
    else:
        surface = QtGui.QSurfaceFormat()
        surface.setRenderableType(QtGui.QSurfaceFormat.OpenGL)
        surface.setProfile(QtGui.QSurfaceFormat.CoreProfile)
        surface.setVersion(4, 5 if (ROOT / "Source/opengl.py").exists() else 6)
        surface.setSamples(4)
        surface.setOption(QtGui.QSurfaceFormat.DebugContext)
    emit("requested", format=fmt(surface))
    probes = [ClearProbe(surface)]
    if WRAPPER or "OpenGL.platform" in sys.modules:
        probes.append(ClearProbe(surface, loader=True))
    # Probes shown AFTER the application's first paint, so the baseline first
    # frame is not affected by our extra contexts.
    if window:
        window.show()

    def show_probes():
        for probe in probes:
            probe.show()

    def capture():
        for probe in probes:
            image = probe.grabFramebuffer()  # Qt resolves multisampling.
            emit("clear-readback", loader="application" if probe.loader else "Qt",
                 null=image.isNull(), rgba=image.pixelColor(image.width() // 2, image.height() // 2).getRgb()
                 if not image.isNull() else None, framesSwapped=probe.frames)
        emit("summary", shaders=COUNTS, widgets=len(WIDGETS), failedWidgets=len(FAILED_WIDGETS),
             note="Readback proves FBO contents, not compositor presentation. Check magenta windows visually.")
    QtCore.QTimer.singleShot(1000, show_probes)
    QtCore.QTimer.singleShot(2500, capture)
    if args.seconds:
        QtCore.QTimer.singleShot(max(3, args.seconds) * 1000, app.quit)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
