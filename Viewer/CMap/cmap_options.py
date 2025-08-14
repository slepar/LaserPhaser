import os

# Define your initial colormap dictionary

custom_cmap_dir = os.path.dirname(os.path.abspath(__file__))

def refresh_cmap_options():
    cmap_dict = {
        "BINARY": "binary",
        "JET": "jet",
        "GREYS": "Greys",
        "VIRIDIS": "viridis",
        "PLASMA": "plasma",
        "INFERNO": "inferno",
        "MAGMA": "magma",
        "CIVIDIS": "cividis",
        "CUSTOM": "custom",
    }
    for file_name in os.listdir(custom_cmap_dir):
        if file_name.endswith(".json"):  
            cmap_name = os.path.splitext(file_name)[0].upper()
            cmap_dict[cmap_name] = file_name
    return cmap_dict

cmap_dict = refresh_cmap_options()