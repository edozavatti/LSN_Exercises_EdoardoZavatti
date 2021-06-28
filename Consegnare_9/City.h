#ifndef __City__
#define __City__

class City {

private:
  double m_x, m_y;
  int m_num;
  
protected:

public:
  // constructors
  City(double, double, int);
  // destructor
  ~City();
  // methods
  double GetX();
  double GetY();
  double Getnum();
  void SetX(double);
  void SetY(double);
  void Setnum(int);
  };

#endif // __Random__
