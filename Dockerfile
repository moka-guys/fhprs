FROM python:3.6
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD ./fh.py /code
ADD ./test/tests_prs.py /code
ADD requirements.txt /code
RUN pip3 install -r requirements.txt
RUN chmod +x /code/fh.py
RUN python -m unittest tests_prs.py
ENTRYPOINT ["python","/code/fh.py"]

