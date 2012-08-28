#include "t_sort.h"

#define QUICK_SORTHRESH 9

template <class T> void swap(T* k, T* j)
{
	T temp;
	temp = *j;
	*j = *k;
	*k = temp;
}

int findpivot(const int& i, const int& j) {return (i+j)/2;}

template <class T> int partition(T* array, int* temp, int l, int r, T pivot)
{
	do {
		while (array[++l] > pivot);
		while (r && array[--r] < pivot);
		swap(array+l, array+r);
		swap(temp+l, temp+r);
	}while (l < r);
	swap(array + l, array + r);
	swap(temp + l, temp + r);
	return l;
}

template <class T> void inssort(T* array, int* temp, int i, int j)
{
	for (int k=i+1; k<=j; k++)
		for (int m=k; (m>i) && (array[m] > array[m-1]); m--)
		{
			swap(array+m, array+m-1);
		    swap(temp+m, temp+m-1);
		}
}

template <class T> void quick_sort(T *array, int* temp, int i, int j)
{
	if(j-i < QUICK_SORTHRESH) inssort(array, temp, i, j);
	else
	{
		int pivotindex = findpivot(i, j);
		swap(array+pivotindex, array+j);
		swap(temp+pivotindex, temp+j);
		int k = partition(array, temp, i-1, j, array[j]);
		swap(array+k, array+j);
		swap(temp+k, temp+j);
		if((k-i)>1) quick_sort(array, temp, i, k-1);
		if((j-k)>1) quick_sort(array, temp, k+1, j);
	}
}

void d_quick_sort(double* array, int* temp, int i, int j)
{
	quick_sort(array, temp, i, j);
}

