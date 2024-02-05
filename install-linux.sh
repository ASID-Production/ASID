# prepare vne
cd VnE
python3 -m venv VisApp
source VisApp\bin\activate
pip install -r requirements.txt
cd .. 


# prepare rsd
cd api_database 
python3 -m venv venv 
venv/Scripts/activate.bat 
pip install -r requirements.txt 
cd django_project 
python3 manage.py migrate 
cd ../.. 

# building dlls
mkdir build-dir
cd build-dir
cmake ..
make

# final cleanup 
rm -rf cpplib 
rm -rf .github 
rm -rf .git 
rm -f CMakeLists.txt 
rm -f CMakeSettings.json 
rm -f .gitignore* 
