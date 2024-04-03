# prepare vne
cd VnE
python3 -m venv VisApp 
source VisApp/bin/activate
pip install -r requirements.txt 
cd .. 

# building dlls
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

# prepare rsd
cd api_database 
python3 -m venv venv 
source venv/bin/activate
pip install -r requirements.txt 
cd ./django_project
python manage.py makemigrations
python manage.py migrate
python manage.py collectstatic
cd ../.. 

# final cleanup 
rm -rf cpplib 
rm -rf .github 
rm -rf .git 
rm -rf build
rm -f CMakeLists.txt 
rm -f CMakeSettings.json 
rm -f .gitignore* 
rm -f VnE/Source/Extensions/ChemPackSource/GenBonds.dll
rm -f VnE/Source/Extensions/ChemPackSource/CMakeLists.txt
rm -f VnE/Source/Extensions/ChemPackSource/dllmain.cpp
rm -f api_database/django_project/modules/c_modules/cpplib.dll
