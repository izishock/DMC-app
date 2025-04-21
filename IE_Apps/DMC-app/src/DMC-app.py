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
global s

# Create semaphore for Data access
Semaphore_Data = Semaphore()
Semaphore_s = Semaphore()

#============================
# MQTT Client
#============================

class mqttclient:

    # Initialize mqttclient object
    def __init__(self,clientname,broker, user, password, metadata_topic, read_topic, write_topic, plc_connection):
        self.IDDict = {}
        self.broker = broker
        self.mqtt_user = user
        self.mqtt_password = password
        self.metadata_topic = metadata_topic
        self.read_topic = read_topic
        self.write_topic = write_topic
        self.plc_connection = plc_connection
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

                Semaphore_Data.acquire()
                Data[name] = msg["vals"][i].get("val", None)
                Semaphore_Data.release() 

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
        self.client.connect(self.broker, 9883, keepalive=60)
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

def read_parameter(jsonfile):
    with open(jsonfile) as params:
        data = json.load(params)
        return data
    
#============================
# Acquire step response sequence
#============================
    
def StepResponseAcquire(my_mqtt_client, write_topic, Sampling_Time, D):
    global Data
    global s
    PLC_WRITE_DATA = PLC_DATA_FORMAT

    CV = np.zeros(2)

    for i in range (30, 101, 10):
        i_idx = int((i-30)/10)

        PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get("Data_DB.Power")
        PLC_WRITE_DATA['vals'][0]['val'] = i
        PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
        my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)

        k_start = 0

        for j in range (70, 81, 10):
            j_idx = int((j-30)/10)
            
            for k in range(k_start, 2*D):
                Semaphore_Data.acquire()
                temp_Data = Data
                Semaphore_Data.release()

                if k == 0 or k == D:
                    alpha = temp_Data["Data_DB.Flow_max"]/(100*60*temp_Data["Data_DB.Volume"])
                    beta = temp_Data["Data_DB.Power_Nominal"]*i/(100*temp_Data["Data_DB.Specific_Heat"]*temp_Data["Data_DB.Density"]*temp_Data["Data_DB.Volume"])
                    G = -0.0002347*i + 1.012

                    PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get("Data_DB.OR")
                    if k == 0:
                        CV[0] = 1/(alpha/(beta*G)*(j-5) - alpha/beta*temp_Data["Data_DB.Temperature_In"])
                        PLC_WRITE_DATA['vals'][0]['val'] = CV[0]
                    elif k == D:
                        CV[1] = 1/(alpha/(beta*G)*(j+5) - alpha/beta*temp_Data["Data_DB.Temperature_In"])
                        PLC_WRITE_DATA['vals'][0]['val'] = CV[1]
                    PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
                    my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)

                if k == D:
                    Temperature_Norm = temp_Data["Data_DB.Temperature_Out"]

                if k >= D:
                    k_idx = k-D

                    if k == int(1.8*D):
                        pass

                    Semaphore_s.acquire()
                    s[j_idx,i_idx,k_idx] = (temp_Data["Data_DB.Temperature_Out"] - Temperature_Norm) / (CV[1] - CV[0])
                    Semaphore_s.release()

                if temp_Data["Data_DB.Manual"] == 0:
                    return 1

                sleep(Sampling_Time)

            CV[0] = CV[1]

            db = sqlite3.connect("./Step_Response.db")
            cursor = db.cursor()

            where_conditions = "T_SP = ? AND Power = ?"
            query = f"SELECT * FROM Step_Response WHERE {where_conditions};"
            cursor.execute(query, tuple([j,i],))
            results = cursor.fetchall()

            if len(results) == 0:
                SR_json = json.dumps(s[j_idx,i_idx].tolist())
                insert_placeholders = ", ".join(["?" for _ in range(3)])
                query_insert = f"INSERT INTO Step_Response (T_SP, Power, SR) VALUES ({insert_placeholders});"
                temp_tuple = list([j,i]) + [SR_json]
                cursor.execute(query_insert, tuple(temp_tuple,))
                db.commit()
            
            k_start = D
    return 0

#============================
# DMC Parameters calculation
#============================

# Calculate ke and Ku parameters
def DMCRegulatorParameters(s, D, N, N1, Nu, l):

    [rows, cols] = s.shape[0], s.shape[1]
    ke = np.zeros((rows, cols))
    ku = np.zeros((rows, cols, D))

    for m in range(cols):
        for n in range(rows):
            # Calculate MP matrix
            MP = np.zeros((N-N1[n][m], D))
            for i in range(0,N-N1[n][m]):
                for j in range(0,D-1):
                    if i+j+N1[n][m] < D:
                        MP[i,j] = s[n][m][i+j+N1[n][m]] - s[n][m][j]
                    else:
                        MP[i,j] = s[n][m][D-1] - s[n][m][j]
                        
            # Calculate M matrix
            M = np.zeros((N-N1[n][m], Nu))
            for i in range(0,Nu):
                for j in range(i,N-N1[n][m]):
                    M[j,i] = s[n][m][j+N1[n][m]-i-1]
                    
            # Calculate ke and Ku
            M_temp = M.T @ M
            eye_size = len(M_temp)

            K_temp = np.linalg.inv(M_temp + l[n][m]*np.eye(eye_size))
            K = K_temp @ M.T

            ku[n,m] = K[0] @ MP
            ke[n,m] = np.sum(K[0])

    return ke, ku

#============================
# DMC Parameters calculation for specific Setpoint
#============================

def CalculateParameters(Temperature_Setpoint, Power, Ke_array, Ku_array):

    Temperature_Setpoint_array = np.array([30, 40, 50, 60, 70, 80])
    Power_array = np.array([30, 40, 50, 60, 70, 80, 90, 100])

    Power_indices = np.where(Power_array >= Power)[0]
    Power_idx = Power_indices[0]
    Temperature_Setpoint_indices = np.where(Temperature_Setpoint_array >= Temperature_Setpoint)[0]
    Temperature_Setpoint_idx = Temperature_Setpoint_indices[0]

    x1 = Power_array[Power_idx - 1]
    x2 = Power_array[Power_idx]
    y1 = Temperature_Setpoint_array[Temperature_Setpoint_idx - 1]
    y2 = Temperature_Setpoint_array[Temperature_Setpoint_idx]

    Q11 = Ke_array[Temperature_Setpoint_idx - 1, Power_idx - 1]
    Q12 = Ke_array[Temperature_Setpoint_idx - 1, Power_idx]
    Q21 = Ke_array[Temperature_Setpoint_idx, Power_idx - 1]
    Q22 = Ke_array[Temperature_Setpoint_idx, Power_idx]
    Ke = 1/((x2-x1)*(y2-y1)) * (Q22*(Power-x1)*(Temperature_Setpoint-y1) + Q21*abs(x2-Power)*abs(Temperature_Setpoint-y1) + Q12*abs(Power-x1)*abs(y2-Temperature_Setpoint) + Q11*abs(x2-Power)*abs(y2-Temperature_Setpoint))

    D = Ku_array.shape[2]
    Ku = np.zeros(D)
    for i in range(D):
        Q11 = Ku_array[Temperature_Setpoint_idx - 1, Power_idx - 1, i]
        Q12 = Ku_array[Temperature_Setpoint_idx - 1, Power_idx, i]
        Q21 = Ku_array[Temperature_Setpoint_idx, Power_idx - 1, i]
        Q22 = Ku_array[Temperature_Setpoint_idx, Power_idx, i]
        Ku[i] = 1/((x2-x1)*(y2-y1)) * (Q22*(Power-x1)*(Temperature_Setpoint-y1) + Q21*(x2-Power)*(Temperature_Setpoint-y1) + Q12*(Power-x1)*(y2-Temperature_Setpoint) + Q11*(x2-Power)*(y2-Temperature_Setpoint))

    return Ke, Ku

#============================
# Main Function
#============================

def main():

    # Establish connection with MQTT broker using parameters from /cfg-data/param.json
    global Data
    global s

    PLC_WRITE_DATA = PLC_DATA_FORMAT

    params = read_parameter("./cfg-data/param.json")

    mqtt_broker_server = params["MQTT_Broker"]
    mqtt_user = params["User"]
    mqtt_password = params["Password"]
    meta_data_topic = params["Metadata_Topic"]
    read_topic = params["Read_Topic"]
    write_topic = params["Write_Topic"]
    variables = params["Variables"]
    variable_control = params["Variable_Control"]
    N = params["DMC_Parameters"][0].get("N")
    Nu = params["DMC_Parameters"][0].get("Nu")
    D = params["DMC_Parameters"][0].get("D")
    N1 = params["DMC_Parameters"][0].get("N1")
    l = params["DMC_Parameters"][0].get("l")
    Sampling_Time = params["DMC_Parameters"][0].get("Sampling_Time")/10

    Data = dict.fromkeys(variables, None)
    s = np.zeros((6, 8, D))

    connection_name = write_topic.split("/")
    connection_name = connection_name[len(connection_name)-1]

    my_mqtt_client = mqttclient('dmc', mqtt_broker_server, mqtt_user, mqtt_password, meta_data_topic,
                                read_topic, write_topic, connection_name)
    my_mqtt_client.start_connection()

    Setpoint_prev = np.zeros(2)

    Semaphore_Data.acquire()
    temp_Data = Data
    Semaphore_Data.release()

    sleep(Sampling_Time)

    while True:
        Setpoint_prev = [temp_Data["Data_DB.Temperature_SetPoint"], temp_Data["Data_DB.Power"]]
        Semaphore_Data.acquire()
        temp_Data = Data
        Semaphore_Data.release()
        Setpoint = [temp_Data["Data_DB.Temperature_SetPoint"], temp_Data["Data_DB.Power"]]

        if temp_Data["Data_DB.Manual"] == 1:
            # error = StepResponseAcquire(my_mqtt_client, write_topic, Sampling_Time, D)
            # print("Error code = " + str(error))

            db = sqlite3.connect("./Step_Response.db")
            cursor = db.cursor()

            for i in range(30, 101, 10):
                for j in range(30, 81, 10):
                    where_conditions = "T_SP = ? AND Power = ?"
                    query = f"SELECT * FROM Step_Response WHERE {where_conditions};"
                    cursor.execute(query, tuple([j,i],))
                    results = cursor.fetchall()

                    i_idx = int((results[0][1]-30)/10)
                    j_idx = int((results[0][0]-30)/10)
                    s[j_idx,i_idx] = json.loads(results[0][2])

            [Ke_array, Ku_array] = DMCRegulatorParameters(s, D, N, N1, Nu, l)

            [Ke, Ku] = CalculateParameters(Setpoint[0], Setpoint[1], Ke_array, Ku_array)
            
            PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get("DMC_Parameters_DB.Ke")
            PLC_WRITE_DATA['vals'][0]['val'] = Ke

            for i in range(0,D):
                index = my_mqtt_client.IDDict.get(f"{"DMC_Parameters_DB.Ku"}[{i}]")
                PLC_WRITE_DATA["vals"].append({"id": index, "val": Ku[i]})
            
            PLC_WRITE_DATA["vals"].append({"id": my_mqtt_client.IDDict.get("Data_DB.Manual"), "val": 0})
            PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
            my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)

        if Setpoint != Setpoint_prev:
            [Ke, Ku] = CalculateParameters(Setpoint[0], Setpoint[1], Ke_array, Ku_array)

            PLC_WRITE_DATA['vals'][0]['id'] = my_mqtt_client.IDDict.get("DMC_Parameters_DB.Ke")
            PLC_WRITE_DATA['vals'][0]['val'] = Ke

            for i in range(0,D-1):
                index = my_mqtt_client.IDDict.get(f"{"DMC_Parameters_DB.Ku"}[{i}]")
                PLC_WRITE_DATA["vals"].append({"id": index, "val": Ku[i]})
            
            PLC_WRITE_DATA["vals"].append({"id": my_mqtt_client.IDDict.get("Data_DB.Manual"), "val": 0})
            PLC_WRITE_DATA_Str = json.dumps(PLC_WRITE_DATA)
            my_mqtt_client.client.publish(write_topic, PLC_WRITE_DATA_Str)

        sleep(Sampling_Time)
    
if __name__ == "__main__":
    main()