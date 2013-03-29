#ifndef _TIMER_H_
#define _TIMER_H_

#include <ctime>
#include <string>

/**
   @class Timer
   @brief A Timer

   This class is a Timer.
 */

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

/**
   @fn Timer::Timer(void)
   Instantiates a new Timer
   **

   @fn Timer::~Timer
   Deletes this Timer
   **

   @fn Timer::start
   Starts this Timer
   **

   @fn Timer::stop
   Stops this Timer

   @note If this Timer hasn't been started,
   an Exception is thrown
   **

   @fn Timer::time
   @return Returns the elapsed time between
   Timer::start() and Timer::stop()

   @note If this Timer hasn't been started
   @em or stoped, an Exception is thrown
   **

   @fn Timer::unit
   @return Returns a string with the @em unit
   used by Timer::time()
 */

#endif
