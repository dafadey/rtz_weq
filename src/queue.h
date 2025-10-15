#pragma once

#include <vector>
#include <mutex>

struct busyQueue {
  enum estate{FREE, BOOKED, WORKING};
  busyQueue();
  
  int thCount = 1;
  
  std::mutex sharedMutex; // this is used to stop client thread. any worker that is ready to get the task unlocks this mutex 
  std::vector<int> slotState; // 0 free, 1 booked, 2 working
  std::mutex slotStateMutex;
  int busyFindSlot();
  bool addTask(int id, void (*f)(void* args, int id), void* args);
  //first you find slot possibly waiting for empty one. the slot will booked after that
  //then YOU prepare data using id of free slot
  //then you run task on this slot. if it is booked it will become busy and 
  bool empty();
};
