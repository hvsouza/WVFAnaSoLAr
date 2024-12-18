#!/usr/bin/env sh
# ________________________________________ #
# Author: Henrique Souza
# Filename: load_ana_class.sh
# Created: 2024
# ________________________________________ #
#
# insert this line to .bashrc
# alias myclass="source __USER_PATH__/load_ana_class.sh"
# to use the sample (analyzer) class, go to any folder which have "analyzed.root" file and run:
# myclass
# If you want to execute a specific file instead, you can run
# myclass {path_to_file}/your_file.root 
# If you want to change old data, you need to give the number of points per waveform as following:
# myclass {file.root} {# of pts}
# If you want just to load the class without any data, you can run
# myclass no
#
# Using `-c` flag will compile the code, which will run much faster
# Add this to your ~/.rootlogon.C:
# gInterpreter->AddIncludePath("__USER_PATH__");
# So scripts are automatically loaded

#!/bin/bash

# Set default values for the parameters
file="analyzed.root"
compile_flag=""

# Process the command-line arguments using getopt
TEMP=$(getopt -o c -- "$@") #
eval set -- "$TEMP"

# Loop through the options and arguments
while true; do
    case "$1" in
        -c)
            compile_flag="+"
            shift # this changed the position 1, 2, 3 to 2, 3 and so on
            ;;
        --) # reach end if parameters
            shift
            ;;
        *)
            if [[ -n "$1" ]]; then
                file="$1"
            fi
            shift
            break
            ;;
    esac
done

# Assign the remaining positional parameters
if [[ $# -ge 1 ]]; then
    echo "Not sure what you want... "
    return 1
fi


try_changing_ext() {
    # local file="$1" # dont need this, but interesting
    new_file="${file%.*}.root"
    if [[ $new_file == '..root' ]]; then #didnt found any
        new_file="$file.root"
    fi

    if [[ -f $new_file ]]; then
        echo "Root file found with the same name, trying this one..."
        file=$new_file
    else
        return 1
    fi
}

check_file_extension() {
    local file="$1"
    local root_ext="${file##*.}"

    if [[ $root_ext != "root" ]]; then
        if ! try_changing_ext; then
            echo "ERROR: The file must have a 'root' extension."
            return 1
        fi
    fi
}

search_for_analyzed_file() {
    local dir="$1"

    if [[ -f "$dir/analyzed.root" ]]; then
        echo "$dir/analyzed.root"
    fi
}

process_file() {
    if [[ $file == "no" ]]; then
       echo "Loading only the class"
       return
    fi
    if [[ -f "$file" ]]; then
        check_file_extension "$file" || return 1
    elif [[ -d "$file" ]]; then
        echo "Searching the directory $file for 'analyzed.root' file..."
        local analyzed_file=$(search_for_analyzed_file "$file")
        if [[ -n "$analyzed_file" ]]; then
            file="$analyzed_file"
        else
            echo "ERROR: no 'analyzed.root' file was found, please provide the file name"
            return 1
        fi
    else
        if ! try_changing_ext; then
            echo "Invalid input: $file"
            return 1
        fi
    fi
    echo "File: $file"
}

# File/Path validadtion:
process_file || return 1


if [ $file == 'no' ]; then
    eval 'root -e ".L __USER_PATH__/WVFAnaSoLAr.h$compile_flag"'
else
    eval 'root -e ".L __USER_PATH__/WVFAnaSoLAr.h$compile_flag" -e "ANALYZER s(\"s\")" -e "s.setAnalyzer(\"$file\")"'
fi


