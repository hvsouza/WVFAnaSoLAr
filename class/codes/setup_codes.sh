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

cp codes/.templates/template_load_my_class.sh load_my_class.sh

sed -i "s|__USER_PATH__|$class_dir|g" load_my_class.sh
