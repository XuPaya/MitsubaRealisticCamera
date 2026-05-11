import time
import mitsuba as mi
import numpy as np

mi.set_variant('cuda_ad_spectral')

import RealisticCamera

start = time.time()
# scene = mi.load_file(
#     'scenes/kitchen/scene_realistic_fisheye.xml'
# )
scene = mi.load_file(
    'scenes/debug/scene.xml'
)
image = mi.render(scene, spp=512)

out = 'outputs/debug_scene.exr'
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