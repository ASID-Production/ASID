
rem prepare sca
cd VnE
call python -m venv VisApp
call VisApp\Scripts\activate.bat
pip install -r requirements.txt
cd ..


rem prepare rsd
cd api_database
call python -m venv venv
call venv\Scripts\activate.bat
pip install -r requirements.txt
pip install psycopg2
cd django_project
call python manage.py makemigrations
call python manage.py migrate
call python manage.py collectstatic
cd ..\..


rem final cleanup
rmdir /S /Q cpplib
rmdir /S /Q docs
rmdir /S /Q .github
rmdir /S /Q .git
del /F /Q README.md
del /F /Q CMakeLists.txt
del /F /Q CMakeSettings.json
del /f /Q .gitignore*
