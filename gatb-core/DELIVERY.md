--------------------------------------------------------------------------------
# HOW TO DELIVER A NEW RELEASE ?

A delivery is made of two files:
* binary archive
* source archive

It is possible to use CMake for building such a delivery. When you invoke 'cmake', 
you have to specify the version number like this:
    
    cmake -DMAJOR=3 -DMINOR=4 -DPATCH=11 ..

Then, you can use a 'delivery' target with:

    make delivery

This target will both produce the binary and source archives and then upload them on 
the GForge server.

IMPORTANT: you have to check by yourself if the version doesn't already exist on 
the server (see below to know what are the existing releases). If you provide an
already existing version number, the previous release on the server will be replaced
by the new one.

You can have the list of available targets for delivery with:

    make delivery_help

--------------------------------------------------------------------------------
# HOW TO KNOW THE EXISTING RELEASES ?

There is a target for this:
    make delivery_dump
I