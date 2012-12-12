#ifndef _TIMER_H_
#define _TIMER_H_

#include <ctime>
#include <string>

class Timer{
 private:
  struct timespec start_s, stop_s;
  double elapsedTime;
  
  bool isLaunched;
  bool hasBeenStoped;

 public:
   Timer(void);
  ~Timer(void);
  
  void        start(void);
  void        stop(void);
  double      time(void);
  std::string unit(void);

 private:
  double timeDiff(void);
};

#endif
