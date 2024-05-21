from OpenGL.GL import *
import ctypes


vars = {}


def get_f(val):
    return val


def get(name):
    func = vars.get(name, None)
    if func is None:
        return 'No such var'
    else:
        return func()


def add_watch(name, val):
    vars[name] = lambda: get_f(val)
    return name


def add_watch_with_ref(ref, name):
    func = lambda: ref.__getattr__(name)
    name_var = 'var_with_ref_' + str(name)
    vars[name_var] = func
    return name_var


def track_buffer_data(buffer, bytes):
    float32_data = (ctypes.c_float * bytes//4)()
    void_ptr = ctypes.c_void_p(ctypes.addressof(float32_data))
    glGetNamedBufferSubData(buffer, 0, bytes, void_ptr)
    ret = [float32_data[x] for x in range(bytes//4)]
    return ret


def track_shader_data(shader_data):
    ret = []
    for vbo in shader_data.VBOs:
        float32_data = (ctypes.c_float * vbo['size'])()
        void_ptr = ctypes.c_void_p(ctypes.addressof(float32_data))
        glGetNamedBufferSubData(vbo['id'], 0, vbo['size_bytes'], void_ptr)
        data = [float32_data[x] for x in range(vbo['size'])]
        ret.append(data)
    return ret


def add_buffer_data(buffer, bytes):
    name = 'buffer_data_' + str(buffer)
    vars[name] = lambda: track_buffer_data(buffer, bytes)
    return name


def add_shader_data(shader_data):
    name = 'shader_data_' + str(shader_data)
    vars[name] = lambda: track_shader_data(shader_data)
    return name
