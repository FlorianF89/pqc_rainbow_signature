load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

cc_library(
    name = "timing-functions",
    srcs = ["timing-functions.c"],
    hdrs = ["timing-functions.h"],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "signature",
    srcs = ["signature_system.c"],
    deps = [
        ":timing-functions",
        "//src/lib:rainbow",
    ],
)
