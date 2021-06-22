## Path Tracing Engine
Path tracing. In this engine, you can make small scenes from simple objects, choose a good angle, wait for the image to be processed and then save the result. You can add more objects, spheres, boxes, ellipsoids, and hexagonal prisms. All you need to do is open the shader.frag in any editor and go down to find a list of objects, there are comments that will help you understand what is responsible for what. And then go up and increase the number of objects. To adjust the quality, change the parameters at the beginning of the shader.frag file. The engine has 2 types of anti-aliasing to smooth out the corners of objects and get more beautiful frames: depth of field, additional anti-aliasing.

NOTE: On some GPUs this may not work. It will be fixed in next updates.
