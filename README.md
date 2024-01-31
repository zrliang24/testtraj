# testtraj

This package requires libgfortran.so.3 and g++-6. Install using *bionic universe* on Ubuntu 20.04, and *not* the focal universe.

If experiencing issues with the get_traj() function, check that the following are true:

-   The hyts_std executable file is in your ExecDir. If you did not specify an ExecDir, make sure that the hyts_std file is in your working directory.

-   The hyts_std executable file has permissions to execute as program. To verify this, go to the properties of the hyts_std file and check the permissions. Ensure that "Allow executing file as program" is checked.
