#pragma once

#ifdef _MSC_VER
#define DCW_BEGIN           \
    __pragma(warning(push)) \
        __pragma(warning(disable : 4305))
#define DCW_END __pragma(warning(pop))
#else
#ifdef __GNUC__ || __clang__
#define DCW_BEGIN               \
    _Pragma("GCC warning push") \
        _Pragma("GCC diagnostic ignored \"-Wconversion\"")
#define DCW_END _Pragma("GCC warning pop")
#endif
#endif

#ifdef DCW_BEGIN
#define DISABLE_CONVERSION_WARNING_BEGIN DCW_BEGIN
#define DISABLE_CONVERSION_WARNING_END DCW_END
#else
#define DISABLE_CONVERSION_WARNING_BEGIN
#define DISABLE_CONVERSION_WARNING_END
#endif
