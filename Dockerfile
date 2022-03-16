FROM python:3.7
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD ./fh.py /code
ADD requirements.txt /code
RUN pip3 install -r requirements.txt
RUN chmod +x /code/fh.py
ENTRYPOINT ["python","/code/fh.py"]

