from matplotlib.colors import LinearSegmentedColormap

cyan_magenta = LinearSegmentedColormap.from_list(
    "cyan_magenta",
    [
        "#00FFFF",
        "#FFFFFF",
        "#FF00FF"
    ],
    N=256
)

cmap_cyan_orange = LinearSegmentedColormap.from_list( # 2nd
    "cyan_white_orange",
    [
        "#00BFFF",   # negative
        "#FFFFFF",   # zero
        "#FF8C00"    # positive
    ],
    N=256
)

cmap_teal_pink = LinearSegmentedColormap.from_list(
    "teal_white_pink",
    [
        "#008C95",
        "#FFFFFF",
        "#E83E8C"
    ],
    N=256
)

cmap_black_red = LinearSegmentedColormap.from_list( # 3rd
    "black_white_red",
    [
        "#000000",
        "#FFFFFF",
        "#D62728"
    ],
    N=256
)