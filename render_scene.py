import time
import mitsuba as mi
import drjit as dr
import numpy as np
from tqdm import tqdm

mi.set_variant('cuda_ad_spectral')

import RealisticCamera

start = time.time()
scene = mi.load_file(
    'scenes/kitchen/scene_realistic_fisheye.xml'
)
# scene = mi.load_file(
#     'scenes/debug/scene.xml'
# )

n_trials = 16

image = None
dr.disable_grad() 
for i in tqdm(range(n_trials)):
    image_i = mi.render(scene, spp=128, seed=i)
    if image is None:
        image = image_i
    else:
        image += image_i
image /= n_trials
out = 'outputs/kitchen_spectral.exr'
mi.Bitmap(image).write(out)

arr = np.array(image)
print('wrote', out)
print(
    'shape', arr.shape,
    'mean', float(arr.mean()),
    'max', float(arr.max()),
    'nonzero', int(np.count_nonzero(arr)),
    'seconds', time.time() - start,
)