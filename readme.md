Old Teaching Code (MATLAB)
==============

This consists of two implementations of the BLP (1995) demand estimation algorithm that I use for teaching purposes only.

**This is NOT production code**

If you want to estimate demand on a large-scale problem I suggest you look at my other repositories.

There are three main files.

blpsimple.m
--------------
This is designed to be a clear and compact version of a simple BLP code for teaching purposes. 
The focus is on expositional clarity and it does not include several "tricks" to speed things up. 
It only handles the most simple cases of independent normal random coefficients. It does not handle other cases.

blpparallel.m
--------------
This is designed to be a clear and compact version of a simple BLP code for teaching purposes. 
The focus is on expositional clarity and it does not include several "tricks" to speed things up. 
It only handles the most simple cases of independent normal random coefficients. It does not handle other cases.

It does show how to parallelize the standard nested fixed point algorithm.

generatedata.m
--------------
This generates some sample data for the above code so that the user understands the format of the data.
