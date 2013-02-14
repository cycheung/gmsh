#ifndef _ELEMENTDATA_H_
#define _ELEMENTDATA_H_

#include "GroupOfDof.h"

class ElementData{
 private:
  bool hasGod;

  const GroupOfDof* god;

 public:
   ElementData(void);
  ~ElementData(void);

  void setGroupOfDof(const GroupOfDof& god);

  const GroupOfDof& getGroupOfDof(void) const;
};

#endif
