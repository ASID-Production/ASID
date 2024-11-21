from OpenGL.GL import *
import ctypes
import time
import sys


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


class Timer:

    class Tim:
        def __init__(self):
            self.tot_time = 0
            self.counting = False
            self.start_tim = 0

        def start(self):
            if not self.counting:
                self.counting = True
                self.start_tim = time.perf_counter()

        def stop(self):
            if self.counting:
                self.tot_time += time.perf_counter() - self.start_tim
                self.counting = False
            else:
                return

    def __init__(self):
        self.groups = {'default': {}}
        self.method_class = {}
        self.class_timers = {}

    def getGroup(self, group):
        if group is None:
            dict = self.groups['default']
        else:
            dict = self.groups.get(group, None)
            if dict is None:
                dict = {}
                self.groups[group] = dict
        return dict

    def time_func(self, func, group=None):
        dict = self.getGroup(group)
        dict[func] = self.Tim()

        vals = func.__qualname__.split('.')[:-1]
        vals = '.'.join(vals)
        self.method_class[func] = vals
        def ret(*args, **kwargs):
            cls_tim = self.class_timers.get(vals, None)
            z = False
            if cls_tim is not None and not cls_tim.counting:
                cls_tim.start()
                z = True
            if dict[func].counting:
                ret_val = func(*args, **kwargs)
            else:
                dict[func].start()
                ret_val = func(*args, **kwargs)
                dict[func].stop()
            if z:
                cls_tim.stop()
            return ret_val
        return ret

    def time_class(self, cls):
        self.class_timers[cls.__name__] = self.Tim()
        return cls

    def print_res(self, group=None):
        if group is None:
            group = 'default'
        dict = self.getGroup(group)
        print('-'*25 + group + '-'*25)
        for item in self.class_timers.items():
            print(f'{item[0]}: {item[1].tot_time}')
        print('*' * 50)
        for item in dict.items():
            print(f'{item[0]}: {item[1].tot_time}')
        print('-'*50)

    def clear(self, group):
        dict = self.getGroup(group)
        for key in dict:
            dict[key] = self.Tim()
        for key in self.class_timers:
            self.class_timers[key] = self.Tim()

TIMER = Timer()