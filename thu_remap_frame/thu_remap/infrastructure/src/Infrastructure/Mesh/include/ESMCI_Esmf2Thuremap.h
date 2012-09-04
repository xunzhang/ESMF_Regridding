/**
 * ESMCI_Esmf2Thuremap.h 
 * Transfer coords related in esmf(after search) into general format
 **/

#ifndef ESMCI_ESMF2THUREMAP_H
#define ESMCI_ESMF2THUREMAP_H

#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>


namespace ESMCI {

struct Point2D {
  double lon;
  double lat;
};

typedef struct Point2D point;

bool operator < (const point &p1, const point &p2) {
  if(p1.lon < p2.lon)
    return true;
  else
    return false;
}

static int compare(const point &a, const point &b) {
  return a.lon < b.lon || a.lat < b.lat; 
}

static int compare_v(const double &a, const double &b) {
  return a < b;
}

static int compare_lon(const point &a, const point &b) {
  return a.lon < b.lon;
}

static int compare_lat(const point &a, const point &b) {
  return a.lat < b.lat;
}


/*
int findBBoxVertexRect(const int &npnts, const double* coords_lon, const double* coords_lat, 
		   double &leftdown_lon, double &leftdown_lat,
		   double &upright_lon, double &upright_lat) {
  
  std::vector<point> points;
  point tmp;
  for(int i = 0; i < npnts; ++i) {
    tmp.lon = coords_lon[i];
    tmp.lat = coords_lat[i];
    points.push_back(tmp);
  }
  
  std::stable_sort(points.begin(), points.end(), compare);
  
  leftdown_lon = points[0].lon;
  leftdown_lat = points[0].lat;
  upright_lon = points[npnts - 1].lon;
  upright_lat = points[npnts - 1].lat;
  
  return 1;
}
*/

class Patch {

public:

  Patch(const int &npnts, const double* coords_lon, const double* coords_lat, const bool* coords_mask) : size(npnts) {
    point tmp;
    for(int i = 0; i < npnts; ++i) {
      tmp.lon = coords_lon[i];
      tmp.lat = coords_lat[i];
      points_set.push_back(tmp);
      points_mask.insert(std::map<point, bool>::value_type(tmp, coords_mask[i]));
	  if (i % 1000 == 0)
		  printf("Patch()\n");
    }
  }
  
  virtual int get_lon_size() {}

  virtual int get_lat_size() {}

  virtual void findVertexScale(const int &size, 
			       double &leftdown_lon, double &leftdown_lat, 
			       double &upright_lon, double &upright_lat) {} 

  virtual void reorder(double* r_coords_lon, double* r_coords_lat) {}

  virtual void get_reorder_mask(bool* mask) {}

protected:
  int size;
  std::vector<point> points_set;
  std::map<point, bool> points_mask;
};

/*
class PatchTotal : public Patch {
public:
private:
};
*/

class PatchRect : public Patch {

public: 

  PatchRect(const int &npnts, const double* coords_lon, const double* coords_lat, const bool* coords_mask) : Patch(npnts, coords_lon, coords_lat, coords_mask) { 
	  compact_by_sort(); 
	 
	  if(v_lon.size() * v_lat.size() != size) {
		std::cout << "v_lon.size is ()" << v_lon.size() << std::endl;
		std::cout << "v_lat.size is ()" << v_lat.size() << std::endl;
		std::cout << "size is ()" << size << std::endl;
	    std::cout << "is not rect, are you kidding me??\n";
	  }
	  
  }
  
  int get_lon_size() { return v_lon.size(); }

  int get_lat_size() { return v_lat.size(); }
    
  bool checkRect(const std::vector<point> &points, const int &size) {
    
    std::map<double, int> lon_num;
    std::map<double, int> lat_num;
    int check_lon = 0, check_lat = 0;
    for(int i = 0; i < size; ++i) {
      
      if(lon_num.count(points[i].lon) == 0) {
        lon_num.insert(std::map<double, int>::value_type(points[i].lon, 1));
        check_lon ++;
      }
      
      if(lat_num.count(points[i].lat) == 0) {
        lat_num.insert(std::map<double, int>::value_type(points[i].lat, 1));
        check_lat ++;
      }
    
    } 
   
    if(check_lon * check_lat == size)  return true;    
    return false;
  
  }

  int findBBoxVertexRect(const int &size, const std::vector<point> &points, 
                         double &leftdown_lon, double &leftdown_lat,
			 double &upright_lon, double &upright_lat) {
  
    std::stable_sort(points_set.begin(), points_set.end(), compare);
  
    leftdown_lon = points_set[0].lon;
    leftdown_lat = points_set[0].lat;
    upright_lon = points_set[size - 1].lon;
    upright_lat = points_set[size - 1].lat;
  
    return 1;
  }
  
  void findVertexScale(const int &size, 
		       double &leftdown_lon, double &leftdown_lat,
		       double &upright_lon, double &upright_lat) { 
    
    if(checkRect(points_set, size) == true)
      findBBoxVertexRect(size, points_set, leftdown_lon, leftdown_lat, upright_lon, upright_lat);
    else {
      std::cout << "It's so sorry your define a wrong type of object." << std::endl;
      exit(1);
    }
  
  }

  void compact_by_sort() {
   
    std::stable_sort(points_set.begin(), points_set.end(), compare_lon); 
    v_lon.push_back(points_set[0].lon);
    
    for(int i = 0; i < size - 1; ++i)
      if(fabs(points_set[i].lon - points_set[i+1].lon) > 1e-10) 
        v_lon.push_back(points_set[i + 1].lon);
    
    std::stable_sort(points_set.begin(), points_set.end(), compare_lat); 
    v_lat.push_back(points_set[0].lat);
      
    for(int i = 0; i < size - 1; ++i) 
      if(fabs(points_set[i].lat - points_set[i+1].lat) > 1e-10) 
        v_lat.push_back(points_set[i + 1].lat);
	std::cout << v_lon.size() << " manimani " << v_lat.size() << std::endl;
    	
  }

  void compact(){

    std::map<double, int> lon_num;
    std::map<double, int> lat_num;

    for(int i = 1; i < size; ++i) {
      if(lon_num.count(points_set[i].lon) == 0) {
        lon_num.insert(std::map<double, int>::value_type(points_set[i].lon, 1));
        v_lon.push_back(points_set[i].lon);
      }

      if(lat_num.count(points_set[i].lat) == 0) {
        lat_num.insert(std::map<double, int>::value_type(points_set[i].lat, 1));
        v_lat.push_back(points_set[i].lat);
      }
    }  
    
    std::sort(v_lon.begin(), v_lon.end(), compare_v);
    std::sort(v_lat.begin(), v_lat.end(), compare_v);
  }  

  void reorder(double* r_coords_lon, double* r_coords_lat) {
   
    for(int i = 0; i < v_lon.size(); ++i) {
      for(int j = 0; j < v_lat.size(); ++j) {
        r_coords_lon[j*v_lon.size()+i] = v_lon[i];
        r_coords_lat[j*v_lon.size()+i] = v_lat[j];
      }
    }
  }

  void get_reorder_mask(bool* mask) {
    
	printf("enter get_reorder_mask.\n");
    point tmp;
    for(int i = 0; i < v_lon.size(); ++i) {
      for(int j = 0; j < v_lat.size(); ++j) {
        tmp.lon = v_lon[i];
		tmp.lat = v_lat[j];
		mask[j * v_lon.size() + i] = points_mask[tmp];
      }
    }
	printf("exit get_reorder_mask.\n");
  }

private:
  std::vector<double> v_lon;
  std::vector<double> v_lat;

}; // end of class PatchRect

} // end of namespace ESMCI


#endif // end of ESMCI_ESMF2THUREMAP_H
