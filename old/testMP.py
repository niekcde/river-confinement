import time
from datetime import datetime as dt
from multiprocessing import Pool
from random import random
from time import sleep


print("project printed immediately.")

def  subMain(x):
    print(x)

def main(x):
    print(x)
    subP = x *[x]
    # time.sleep(x)
    with Pool(2) as p:
        A = p.imap(subMain, subP)
        p.close()
        # # wait for all issued task to complete
        p.join()
    return (x*x, x)


tm1 = dt.now()

if __name__ == '__main__':
    with Pool(2) as p:
        A = p.imap(main, [1,5, 3, 8, 10, 15,3])
        p.close()
        # # wait for all issued task to complete
        p.join()

# print('A: ', A)
# print('P: ', p)
# for i in A:
#     print(i[0], i[1])
# print(dt.now() - tm1)

# SuperFastPython.com
# example of parallel imap() with the process pool


 
# # task executed in a worker process
# def task(identifier):
#     # generate a value
#     value = identifier
#     # report a message
#     print(f'Task {identifier} executing with {value}', flush=True)
#     # block for a moment
#     sleep(value)
#     # return the generated value
#     return value
 
# tm1 = dt.now()
# # protect the entry point
# if __name__ == '__main__':
#     # create and configure the process pool
#     with Pool(2) as pool:
#         # execute tasks in order
#         for result in pool.imap(task, [1,10,2,1,3,10,20,10,5,10, 30,10,2,10,3,10,4,4,5,10]):
#             print(f'Got result: {result}', flush=True)
#     # process pool is closed automatically

# print(dt.now() - tm1)

#%%
import logging
from multiprocessing import Pool

# Configure logging for all processes
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [PID %(process)d - %(processName)s] %(message)s",
)

def nested_worker(task_id):
    """
    Function to be executed by subprocesses of a worker.
    """
    logging.info(f"Nested worker for task {task_id} is running. PID: {os.getpid()}")
    # Simulate some work
    return f"Nested task {task_id} completed by PID {os.getpid()}"

def worker_function(main_task_id):
    """
    Main worker function, which itself spawns subprocesses.
    """
    logging.info(f"Worker {main_task_id} is starting. PID: {os.getpid()}")

    # Simulate a list of subtasks for this worker
    subtasks = [f"{main_task_id}-{i}" for i in range(3)]

    # Use multiprocessing to process these subtasks
    with Pool(processes=2) as pool:
        results = pool.map(nested_worker, subtasks)

    logging.info(f"Worker {main_task_id} completed. Results: {results}")
    return f"Worker {main_task_id} finished"

def main():
    """
    Main process that spawns workers.
    """
    logging.info("Main process is starting.")

    # List of main tasks
    main_tasks = range(4)  # 4 main tasks

    # Use multiprocessing to handle main tasks
    with Pool(processes=2) as pool:
        results = pool.map(worker_function, main_tasks)

    logging.info(f"All tasks completed. Final results: {results}")

if __name__ == "__main__":
    main()
# %%
