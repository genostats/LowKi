#ifndef findPos
#define findPos

template<typename dbvec>
int find_pos(double x, dbvec X) {
  int i = 0;
  int j = X.size() - 1;
  if(x < X[0])
    return -1;
  if(x > X[j])
    return j;
  while(i < j - 1) {
    int k = (i+j)/2;
    if(x < X[k]) 
      j = k; 
    else 
      i = k;
  }
  return i;
}

#endif