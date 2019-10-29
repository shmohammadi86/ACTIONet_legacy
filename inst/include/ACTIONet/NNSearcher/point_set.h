#ifndef POINT_SET_H
#define POINT_SET_H

#include <algorithm>

template <class METRIC> class point_set_base {
protected:
  const long N; // number of points
  point_set_base<METRIC>(const long n) : N(n){};
public:
  ~point_set_base<METRIC>(){};
  inline long size() const { return N; };
  typedef METRIC Metric;
};

template <class METRIC>
class rm_point_set : public point_set_base<METRIC> {

protected:
  const long D; // dimension
  float* matrix_ptr; // points are stored row-major in a C style array
  const METRIC Distance; // a function object that calculates distances
public:
  rm_point_set() = delete;
  rm_point_set(const rm_point_set& from) = delete;
  rm_point_set(const mat& m)
    : point_set_base<METRIC>(m.n_cols), D(m.n_rows), matrix_ptr(new float[m.n_elem]), Distance(){
		for(register int i = 0; i < m.n_elem; i++) 
			matrix_ptr[i] = m(i);
		//memcpy(matrix_ptr, m.memptr(), m.n_elem * sizeof(float));
    };
  // Move constructor here.
  rm_point_set(rm_point_set&& from)
    : point_set_base<METRIC>(from.N), D(from.D), Distance(){
      matrix_ptr = from.matrix_ptr;
      from.matrix_ptr = nullptr;
    };
  ~rm_point_set(){
    delete[] matrix_ptr;
  };
  inline long dimension() const { return D; };

  typedef const float* row_iterator; // pointer that iterates over the elements
  // of one point in the rm_point_set (points are row vectors)

  row_iterator point_begin(const long n) const { return matrix_ptr + n*D; }
  row_iterator point_end(const long n) const {
    return matrix_ptr + (n + 1) * D; // past-the-end
  }

  template <class ForwardIterator>
  inline double distance(const long index1, ForwardIterator vec2) const {
    return Distance(point_begin(index1), point_end(index1), vec2);
  }
  template <class ForwardIterator>
  inline double distance(const long index1, ForwardIterator vec2,
                         const double thresh) const {
    return Distance(point_begin(index1), point_end(index1), vec2, thresh);
  }
  inline double distance(const long index1, const long index2) const {
    return Distance(point_begin(index1), point_end(index1),
                    point_begin(index2));
  }
};

#endif
