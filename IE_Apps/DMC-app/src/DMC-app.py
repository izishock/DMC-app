import paho.mqtt.client as mqtt
import sys
import os
import json
import re
import sqlite3
import numpy as np
from time import sleep
from threading import *
from Create_DB import Create_DB

# Data format for MQTT publishing
PLC_DATA_FORMAT = {"seq": 1, "vals": [{"id": "", "val": ""}]}

# Create global variable for publish/subscribe PLC data
global Data
# Create semaphore for Data access
Semaphore_Data = Semaphore()

#============================
# MQTT Client
#============================

class mqttclient:

    # Initialize mqttclient object
    def __init__(self,clientname,broker, user, password, metadata_topic, read_topic, write_topic, plc_connection, variable_array, number_operating_point):
        self.IDDict = {}
        self.broker = broker
        self.mqtt_user = user
        self.mqtt_password = password
        self.metadata_topic = metadata_topic
        self.read_topic = read_topic
        self.write_topic = write_topic
        self.plc_connection = plc_connection
        self.variable_array = variable_array
        self.number_operating_point = number_operating_point
        self.clientname = clientname
        self.client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2)
        self.client.on_connect = self.on_connect
        self.client.on_subscribe = self.on_subscribe
        self.client.on_message = self.on_message        
        
    # Function to call whenever message is received
    def on_message(self, client, userdata, message):
        global Data

        # Clasify message by topic
        if message.topic == self.metadata_topic:

            # Create dictionary with connections id + tag if message is metadata
            msg = message.payload.decode("utf-8")      
            self.IDDict = getDict(msg, self.plc_connection)
        
        elif message.topic == self.read_topic:

            # Process data if message is standard read data
            msg = json.loads(message.payload.decode("utf-8"))
            msg_length = len(msg["vals"])

            # Iterate through the whole message
            for i in range(msg_length):

                # Find tag name in message
                id = msg["vals"][i].get("id", None)
                name = next((key for key, value in self.IDDict.items() if value == id), None)

                # Put tag value in the correct spot in Data
                for j in range(len(self.variable_array)):
                    if name.startswith(self.variable_array[j]):
                        Semaphore_Data.acquire()
                        Data[j] = int(msg["vals"][i].get("val", None))
                        Semaphore_Data.release()                  
                        break

    # Function to call whenever connection is established
    def on_connect(self, client, userdata, flags, rc, properties):
        print("INFO | Connected to " + self.broker + " MQTT broker")
        print("INFO | Subscribe to " + self.metadata_topic + " topic")
        print("INFO | Subscribe to " + self.read_topic + " topic")
        print("INFO | Subscribe to " + self.write_topic + " topic")
        sys.stdout.flush()
        self.client.subscribe(self.metadata_topic)
        self.client.subscribe(self.read_topic)
        self.client.subscribe(self.write_topic)
      
    # Function to call whenever topic is subscribed
    def on_subscribe(self, client, userdata, mid, granted_qos, properties):
        print("INFO | On subscribe callback executed")
        sys.stdout.flush()

    # Function to call whenever client wants to establish connection
    def start_connection(self):
        print("INFO | Start connection to " + self.broker + " broker")
        sys.stdout.flush()
        self.client.username_pw_set(self.mqtt_user, self.mqtt_password)
        self.client.connect(self.broker, 1883, keepalive=60)
        self.client.loop_start()

#============================
# Help functions
#============================

# Return tag name / tag id dictionary 
def getDict(metadata, connection_name):
    nameIDMap = {}
    json_metadata = json.loads(metadata)

    connections = json_metadata['connections']
    for element in connections:
        datapointobjects = element['dataPoints']

        if (element['name']==connection_name):

            for element in datapointobjects:
                definitions = element['dataPointDefinitions']

                for element in definitions:
                    nameIDMap.update({element['name'] : element['id']})

    return(nameIDMap)

# Check if database exists / if not create the correct one
def check_db_connection(db_path):
    if not os.path.exists(db_path):
        db = Create_DB(db_path)
    else:
        db = sqlite3.connect(db_path)
    return db

# Read json file
def read_parameter(jsonfile):
    with open(jsonfile) as params:
        data = json.load(params)
        return data
    
#============================
# DMC Parameters calculation
#============================

# Calculate ke and Ku parameters
def DMCRegulatorParameters(s, D, N, N1, Nu, l):

    # Calculate MP matrix
    MP = np.zeros((N-N1+1, D-1))
    for i in range(0,N-N1+1):
        for j in range(0,D-1):
            if i+j+N1 < D:
                MP[i,j] = s[i+j+N1] - s[j]
            else:
                MP[i,j] = s[D-1] - s[j]
                
    # Calculate M matrix
    M = np.zeros((N-N1+1, Nu))
    for i in range(0,Nu):
        for j in range(i,N-N1+1):
            M[j,i] = s[j+N1-i-1]
            
    # Calculate ke and Ku
    M_temp = M.T @ M
    eye_size = len(M_temp)

    K_temp = np.linalg.inv(M_temp + l*np.eye(eye_size))
    K = K_temp @ M.T

    ku = K[0] @ MP
    ke = np.sum(K[0])

    return ku, ke

#============================
# Main Function
#============================

def main():
    # Connect to database, establish connection with MQTT broker using parameters from /cfg-data/param.json
    global Data

    PLC_WRITE_DATA = PLC_DATA_FORMAT

    db = check_db_connection("./DMC_Parameters_DB.db")
    cursor = db.cursor()

    params = read_parameter("/cfg-data/param.json")

    mqtt_broker_server = params["MQTT_Broker"]
    mqtt_user = params["User"]
    mqtt_password = params["Password"]
    meta_data_topic = params["Metadata_Topic"]
    read_topic = params["Read_Topic"]
    write_topic = params["Write_Topic"]
    number_operating_point = params["Number_Variables_Operating_Point"]
    variables_operating_point = []
    for i in range(number_operating_point):
        variables_operating_point.append(params["Variables_Operating_Point"][i])
    variable_process_value = params["Variable_Step_Response"]
    variable_control = params["Variable_Control"]
    variable_clear_database = params["Variable_Clear_Database"]
    variable_array = [variable_process_value, variable_clear_database, *variables_operating_point]
    N = params["DMC_Parameters"][0].get("N")
    Nu = params["DMC_Parameters"][0].get("Nu")
    D = params["DMC_Parameters"][0].get("D")
    N1 = params["DMC_Parameters"][0].get("N1")
    l = params["DMC_Parameters"][0].get("l")
    Sampling_Time = params["DMC_Parameters"][0].get("Sampling_Time")

    connection_name = write_topic.split("/")
    connection_name = connection_name[len(connection_name)-1]

    my_mqtt_client = mqttclient('dmc', mqtt_broker_server, mqtt_user, mqtt_password, meta_data_topic,
                                read_topic, write_topic, connection_name, variable_array, number_operating_point)
    my_mqtt_client.start_connection()

    # Create required variables
    Semaphore_Data.acquire()
    Data = np.zeros(2 + number_operating_point)  
    Semaphore_Data.release()

    State = 0
    Ke = 0
    Ku = np.zeros(D)
    Control = 0
    Operating_Variables_Prev = np.zeros(number_operating_point)

    sleep(Sampling_Time)

    Semaphore_Data.acquire()
    Operating_Variables_Prev = Data[2:]
    Semaphore_Data.release()

    # Wait for DMC parameters request
    while True:
        Semaphore_Data.acquire()
        PV = Data[0]
        Clear_Database = Data[1]
        Operating_Variables = Data[2:]
        Semaphore_Data.release()

        # Process clearing database if requested
        if Clear_Database == 1:
            cursor.execute("DELETE FROM Parameters")
            db.commit()
            PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get(variable_clear_database)
            PLC_WRITE_DATA['vals'][0]['val'] = 0
            PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
            my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)
            PLC_WRITE_DATA = {"seq": 1, "vals": [{"id": "", "val": ""}]}

        # Check operating point change
        if State == 0 and (Operating_Variables != Operating_Variables_Prev).any:
            where_conditions = " AND ".join([f"Operating_Point_{i+1} = ?" for i in range(number_operating_point)])
            query = f"SELECT * FROM Parameters WHERE {where_conditions};"
            cursor.execute(query, tuple(Operating_Variables,))
            results = cursor.fetchall()

            if len(results) == 0:
                if Ke == 0:
                    State = 1
                else:
                    State = 2
            else:
                Ke = results[0][number_operating_point]
                Ku = json.loads(results[0][number_operating_point+1])
                State = 0
        Operating_Variables_Prev = Operating_Variables

        if State == 1:
            Control = (0.965 * Operating_Variables[0] - 17) / ((3 * 12000) / (5 * Operating_Variables[1] * 4200 * 1))











        # Check database whether it contains parameters for given operating point
        if temp_Data[0] == 1:
            where_conditions = " AND ".join([f"Operating_Point_{i+1} = ?" for i in range(number_operating_point)])
            query = f"SELECT * FROM Parameters WHERE {where_conditions};"
            cursor.execute(query, tuple(temp_Data[1:number_operating_point+1],))
            results = cursor.fetchall()

            if len(results) == 0:
                PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get(variable_state)
                PLC_WRITE_DATA['vals'][0]['val'] = 2
                PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
                my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)
                PLC_WRITE_DATA = {"seq": 1, "vals": [{"id": "", "val": ""}]}
            else:
                ke_retrieved = results[0][number_operating_point]
                PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get(variable_Ke)
                PLC_WRITE_DATA['vals'][0]['val'] = ke_retrieved

                Ku_retrieved = json.loads(results[0][number_operating_point+1])
                for i in range(0,D-1):
                    index = my_mqtt_client.IDDict.get(f"{variable_Ku}[{i}]")
                    PLC_WRITE_DATA["vals"].append({"id": index, "val": Ku_retrieved[i]})

                PLC_WRITE_DATA["vals"].append({"id": my_mqtt_client.IDDict.get(variable_state), "val": 0})

                PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
                my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)
                PLC_WRITE_DATA = {"seq": 1, "vals": [{"id": "", "val": ""}]}

        # Process parameters calculations using read step response vector and parameters from /cfg-data/param.json
        elif temp_Data[0] == 4:
            s = temp_Data[number_operating_point+1:number_operating_point+1+max(N,D)]

            reg_Param = DMCRegulatorParameters(s, D, N, N1, Nu, l)
            PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get(variable_Ke)
            PLC_WRITE_DATA['vals'][0]['val'] = reg_Param[1]

            for i in range(0,D-1):
                index = my_mqtt_client.IDDict.get(f"{variable_Ku}[{i}]")
                PLC_WRITE_DATA["vals"].append({"id": index, "val": reg_Param[0][i]})
            
            PLC_WRITE_DATA["vals"].append({"id": my_mqtt_client.IDDict.get(variable_state), "val": 0})
            PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
            my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)
            PLC_WRITE_DATA = {"seq": 1, "vals": [{"id": "", "val": ""}]}
            
            # Write parameters to database
            if reg_Param[1] != 0 and len(reg_Param[0]) != 0:
                Ke_json = json.dumps(reg_Param[0].tolist())
                insert_columns = ", ".join([f"Operating_Point_{i+1}" for i in range(number_operating_point)])
                insert_placeholders = ", ".join(["?" for _ in range(number_operating_point+2)])
                query_insert = f"INSERT INTO Parameters ({insert_columns}, Ke, Ku) VALUES ({insert_placeholders});"
                temp_tuple = list(temp_Data[1:1 + number_operating_point]) + [reg_Param[1]] + [Ke_json]
                cursor.execute(query_insert, tuple(temp_tuple,))
                db.commit()

        sleep(1)

if __name__ == "__main__":
    main()