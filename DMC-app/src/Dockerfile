FROM python:3.12

WORKDIR /app

COPY DMC_MQTT.py .
COPY Create_DB.py .

RUN pip install paho-mqtt && \
    pip install numpy

CMD [ "python", "./DMC_MQTT.py"]

