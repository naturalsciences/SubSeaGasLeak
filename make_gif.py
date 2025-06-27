import os
from PIL import Image

# Set the directory containing your PNG files
directory = "tmp/"

# Get all PNG files and sort them alphabetically
images = sorted([f for f in os.listdir(directory) if f.endswith('.png')])
# Open all images
frames = [Image.open(os.path.join(directory, image)) for image in images]

# Save as GIF
frames[0].save('output.gif', format='GIF', append_images=frames[1:], save_all=True, duration=0.25*1000, loop=0)

