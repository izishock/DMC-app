import sqlite3
import json

def read_parameter(jsonfile):
    with open(jsonfile) as params:
        data = json.load(params)
        return data

def Create_DB(db_path):
    db = sqlite3.connect(db_path)
    cursor = db.cursor()

    params = read_parameter("/cfg-data/param.json")

    number_operating_point = params["Number_Variables_Operating_Point"]

    column_names = [f"Operating_Point_{i}" for i in range(1, number_operating_point + 1)]
    column_definitions = ", ".join([f"{col_name} FLOAT" for col_name in column_names])

    table_name = "Parameters"
    create_table_query = f"CREATE TABLE IF NOT EXISTS {table_name} ({column_definitions}, Ke FLOAT, Ku FLOAT);"

    cursor.execute(create_table_query)

    return db