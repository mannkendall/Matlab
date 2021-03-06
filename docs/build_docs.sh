# Custom bash script to generate the mannkendall docs
# Created August 2020; F.P.A. Vogt; frederic.vogt@meteoswiss.ch
#

# Step 0: just clean the existing apidoc rst files and any previous build folder
#rm -rf ./source/modules
#rm -rf ./build
rm -rf ./build/*.html
rm -rf ./build/_sources
rm -rf ./build/_static
rm -rf ./build/.buildinfo

# Step 1: run autodoc to generate all the docstring rst files.
# Force the rewrite of all of them to capture *all* the recent changes.
#sphinx-apidoc -f -M -o ./source/modules/ ../src/

# Delete the superfluous module.rst file that the previous command creates.
#rm -f ./source/modules/modules.rst

# Generate the documentation, storing it in the build directory
sphinx-build -a -b html ./source ./build
