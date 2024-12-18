#!/usr/bin/env sh

current_dir=$( pwd )
cd ../
class_dir=$( pwd )
cd $current_dir

# echo "$class_dir" > userpath.log

files=$(/bin/ls -1 .templates/ | grep .C)

for f in $files; do
    new_file=$(echo $f | sed "s/template_//g")
    # echo "$new_file"
    cp .templates/$f $new_file
    sed -i "s|__USER_PATH__|$class_dir|g" "$new_file"
done


cd ../

cp codes/.templates/template_load_ana_class.sh load_ana_class.sh

sed -i "s|__USER_PATH__|$class_dir|g" load_ana_class.sh

if ! grep -Fq load_ana_class ~/.bashrc
then
    echo "Please, add this to your .bashrc:"
    echo "alias anasol=\"source $class_dir/load_ana_class.sh"
fi

if ! grep -Fq $class_dir ~/.rootlogon.C
then
    echo "\nPlease, add this to your .rootlogon.C:"
    echo "gInterpreter->AddIncludePath(\"$class_dir\");"
fi


