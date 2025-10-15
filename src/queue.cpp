#include "queue.h"

#include <omp.h>
#include <thread>
#include <iostream>

busyQueue::busyQueue() : thCount{/*omp_get_num_threads()*/16}, slotState(thCount, busyQueue::estate::FREE) {
    std::cout << "busyQueue total number of slots is " << thCount << '\n';
}

bool busyQueue::empty()
{
  bool isempty = true;
  for(auto ss : slotState)
  {
    slotStateMutex.lock();
    isempty &= (ss == busyQueue::estate::FREE);
    slotStateMutex.unlock();
  }
  return isempty;
}

int busyQueue::busyFindSlot()
{
  int id=-1;
  while(true) {
    slotStateMutex.lock();
    for(int i=0;i<slotState.size();i++) {
      if(slotState[i] == busyQueue::estate::FREE) {
        id=i;
        break;
      }
    }
    if(id != -1)
      slotState[id] = busyQueue::estate::BOOKED;
    slotStateMutex.unlock();
    if(id != -1)
      break;
    sharedMutex.lock();
  }
  
  return id;
}

bool busyQueue::addTask(int id, void (*f)(void*, int id), void* args)
{
  bool slotIsFree = false;
  slotStateMutex.lock();
  if(slotState[id] == busyQueue::estate::BOOKED) {
    slotState[id] = busyQueue::estate::WORKING;
    slotIsFree = true;
  }
  slotStateMutex.unlock();
  if(!slotIsFree)
    return false;
  
  std::thread th([id, f, args, this]() {
                  (*f)(args, id);
                  slotStateMutex.lock();
                  slotState[id] = busyQueue::estate::FREE;
                  slotStateMutex.unlock();
                  sharedMutex.unlock();
                });
  th.detach();
  return true;
}
