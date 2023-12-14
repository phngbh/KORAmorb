#!/bin/bash

# USAGE: setup_rampdb.sh <mysql_config_file> <your_datbase_name> <local_data_dump>

# Install MySQL
brew install mysql # Make sure you have homebrew installed in your computer

# Start MySQL service
brew services start mysql

# Load MySQL credentials from the configuration file
MYSQL_CONFIG_FILE="$1"

# Replace your_database_name with the desired database name
DATABASE_NAME="$2"

# Create the database
mysql --defaults-extra-file=$MYSQL_CONFIG_FILE -e "CREATE DATABASE $DATABASE_NAME;"

# Import data dump into the created database
DUMP_FILE="$3"  

# Check if the dump file exists
if [ -f "$DUMP_FILE" ]; then
    mysql --defaults-extra-file=$MYSQL_CONFIG_FILE $DATABASE_NAME < $DUMP_FILE
    echo "Data dump imported successfully."
else
    echo "Error: Data dump file not found."
fi
