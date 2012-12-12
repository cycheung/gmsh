#include "Timer.h"
#include "Exception.h"

Timer::Timer(void){
  isLaunched    = false;
  hasBeenStoped = false;
}

Timer::~Timer(void){
}

void Timer::start(void){
  clock_gettime(CLOCK_REALTIME, &start_s); 
  isLaunched = true;
}

void Timer::stop(void){
  // For better accuracy, take Time first ...
  clock_gettime(CLOCK_REALTIME, &stop_s); 
  
  if(!isLaunched)
    throw Exception("Trying to Stop Timer without Starting it");
  
  elapsedTime   = timeDiff();
  hasBeenStoped = true; 
}

double Timer::time(void){  
  if(!hasBeenStoped)  
    throw Exception("Trying to Have Elapsed Time without Stoping Timer");  

  if(!isLaunched)  
    throw Exception("Trying to Have Elapsed Time without Launching Timer");  

  return elapsedTime;
}

std::string Timer::unit(void){
  return std::string("s");
}

double Timer::timeDiff(void){
  if(start_s.tv_nsec < stop_s.tv_nsec)
    return (stop_s.tv_sec - start_s.tv_sec) 
      + ((double)(stop_s.tv_nsec - start_s.tv_nsec) / 1000000000);
  else
    return (stop_s.tv_sec - start_s.tv_sec - 1)
      + ((double)(1000000000 + stop_s.tv_nsec - start_s.tv_nsec) / 1000000000);
}
