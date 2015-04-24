--------------------------------------------------------------------------------
# HOW TO DELIVER A NEW RELEASE ?

A delivery is made of two files:
* binary archive
* source archive

It is possible to use CMake for building such a delivery. When you invoke 'cmake', 
you have to specify the version number like this:
    
    cmake -DMAJOR=1 -DMINOR=2 -DPATCH=11 ..

For better compatibility, you might want to create a static binary, by replacing 
the above command line with:
    
    cmake -Dstatic=1 -DMAJOR=1 -DMINOR=2 -DPATCH=11 ..

Then, you can use a 'delivery' target with:

    make delivery

This target will both produce the binary and source archives and then upload them on 
the GForge server.

You can have the list of available targets for delivery with:

    make delivery_help

--------------------------------------------------------------------------------
# HOW TO KNOW THE EXISTING RELEASES ?

There is a target for this:
    make delivery_dump
I
