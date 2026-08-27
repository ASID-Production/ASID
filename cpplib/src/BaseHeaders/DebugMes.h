#pragma once

#ifdef DEBUG_MESS_ON
#include <concepts>
#include <iostream>
#include <string_view>
#include <utility>

namespace debug {
    struct LogDepth {
        static inline thread_local int level = 0; // thread_local for multithreadting
        static void print_indent() {
            for (int i = 0; i < level; ++i) std::cerr << "  ";
        }
    };

    struct InterfaceGuard {
        std::string_view name;

        explicit InterfaceGuard(std::string_view n) : name(n) {
            LogDepth::print_indent();
            std::cerr << "[START] " << name << std::endl;
            LogDepth::level++;
        }

        ~InterfaceGuard() {
            LogDepth::level--;
            LogDepth::print_indent();
            std::cerr << "[END]   " << name << std::endl;
        }

        InterfaceGuard(const InterfaceGuard&) = delete;
        InterfaceGuard& operator=(const InterfaceGuard&) = delete;
    };

    template <typename F, typename... Args>
        requires std::invocable<F, Args...>
    decltype(auto) execute(std::string_view name, F&& func, Args&&... args) {
        InterfaceGuard scope(name);
        return std::forward<F>(func)(std::forward<Args>(args)...);
    }

    template <typename Obj, typename Method, typename... Args>
    decltype(auto) execute_method(std::string_view method_name, Obj& obj, Method method, Args&&... args) {
        InterfaceGuard scope(method_name);
        return (obj.*method)(std::forward<Args>(args)...);
    }
}

#define DEBUG_CONCAT_HIDDEN(a, b) a##b
#define DEBUG_CONCAT(a, b) DEBUG_CONCAT_HIDDEN(a, b)

#define WITH_LOG(func, ...) debug::execute(#func, func, ##__VA_ARGS__)
#define LOG_INTERFACE_GUARD(name) debug::InterfaceGuard DEBUG_CONCAT(guard_, __LINE__)##(name)
#define WITH_LOG_M(obj, method, ...) \
    debug::execute_method(#method, obj, &std::remove_pointer_t<decltype(&obj)>::method, ##__VA_ARGS__)

#else
#define WITH_LOG(func, ...) func(__VA_ARGS__)
#define WITH_LOG_M(obj, method, ...) (obj.method(__VA_ARGS__))
#define LOG_INTERFACE_GUARD(name) [[maybe_unused]] int dummy_##__LINE__ = 0
#endif