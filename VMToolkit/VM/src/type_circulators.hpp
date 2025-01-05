#ifndef __TYPE_CIRCULATORS_HPP__
#define __TYPE_CIRCULATORS_HPP__


#include "types.hpp"

namespace VMTutorial
{
    // Forward declarations
	
	class HalfEdge;
	
	class Vertex;
	
	class Edge;
	
	class Face;
	
	class Mesh;

	// Vertex circulator
	class VertexCirculator
	{

	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = HalfEdge;
		using difference_type = std::ptrdiff_t;
		using pointer = HalfEdge *;
		using reference = HalfEdge &;

		VertexCirculator() : _start{},
							 _current{},
							 _isEnd{true}
		{
		}

		explicit VertexCirculator(HEHandle he) : _start{he},
														   _current{he},
														   _isEnd{false}
		{
		}
		explicit VertexCirculator(Vertex v) : _start{v.he()},
												        _current{v.he()},
												        _isEnd{false}

		{
		}

		VertexCirculator &operator++()
		{
			_current = _current->pair()->next();
			if (_current == _start)
			{
				_isEnd = true; // Completed full circle
			}
			return *this;
		}

		VertexCirculator operator++(int)
		{
			VertexCirculator temp = *this;
			++(*this);
            // *temp -> *_current -> (HalfEdge &)*HEHandle
			// return *temp;
			return temp;
		}

		reference operator*()
		{
			return *_current;
		}

		pointer operator->()
		{
			return &(*_current);
		}

		bool operator==(const VertexCirculator &other) const
		{
			return (_isEnd && other._isEnd) || (_current == other._current);
		}

		bool operator!=(const VertexCirculator &other) const
		{
			return !(*this == other);
		}

		// Python iterator protocol methods
		VertexCirculator& __iter__() 
		{
			return *this;
		}

		const reference __next__() 
		{
			if (_isEnd) 
			{
				throw py::stop_iteration();
			}
			HEHandle result = _current;
			this->operator++();
			return *result;
		}

		// Begin and end methods for range-based for loop
		VertexCirculator begin() { return *this; }
		VertexCirculator end() { return VertexCirculator(); }

	private:
		HEHandle _start;
		HEHandle _current;
		bool _isEnd;
	};

	// constant Vertex circulator
	class VertexCCirculator
	{

	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = HalfEdge;
		using difference_type = std::ptrdiff_t;
		using pointer = const HalfEdge *;
		using reference = const HalfEdge &;

		VertexCCirculator() : _start{},
							  _current{},
							  _isEnd{true}
		{
		}

		explicit VertexCCirculator(HECHandle he) : _start{he},
															 _current{he},
															 _isEnd{false}
		{
		}

		VertexCCirculator &operator++()
		{
			_current = _current->pair()->next();
			if (_current == _start)
			{
				_isEnd = true; // Completed full circle
			}
			return *this;
		}

		VertexCCirculator operator++(int)
		{
			VertexCCirculator temp = *this;
			++(*this);
			// return *temp;
			return temp;
		}

		reference operator*() const
		{
			return *_current;
		}

		pointer operator->()
		{
			return &(*_current);
		}

		bool operator==(const VertexCCirculator &other) const
		{
			return (_isEnd && other._isEnd) || (_current == other._current);
		}

		bool operator!=(const VertexCCirculator &other) const
		{
			return !(*this == other);
		}
		

		// Begin and end methods for range-based for loop
		VertexCCirculator begin() { return *this; }
		VertexCCirculator end() { return VertexCCirculator(); }

	private:
		HECHandle _start;
		HECHandle _current;
		bool _isEnd;
	};

	// Face Circulator
	class FaceCirculator
	{

	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = HalfEdge;
		using difference_type = std::ptrdiff_t;
		using pointer = const HalfEdge *;
		using reference = const HalfEdge &;

		FaceCirculator() : _start{},
						   _current{},
						   _isEnd{true}

		{
		}
		explicit FaceCirculator(HEHandle he) : _start{he},
														 _current{he},
														 _isEnd{false}

		{
		}
		explicit FaceCirculator(Face f) : _start{f.he()},
												    _current{f.he()},
													_isEnd{false}

		{
		}

		FaceCirculator &operator++()
		{
			_current = _current->next();
			if (_current == _start)
			{
				_isEnd = true; // Completed full circle
			}
			return *this;
		}

		FaceCirculator operator++(int)
		{
			FaceCirculator temp = *this;
			++(*this);
			// return *temp;
			return temp;
		}

		reference operator*()
		{
			return *_current;
		}

		pointer operator->()
		{
			return &(*_current);
		}

		bool operator==(const FaceCirculator &other) const
		{
			return (_isEnd && other._isEnd) || (_current == other._current);
		}

		bool operator!=(const FaceCirculator &other) const
		{
			return !(*this == other);
		}

		// Python iterator protocol methods
		FaceCirculator& __iter__() 
		{
			return *this;
		}

		reference __next__() 
		{
			if (_isEnd) 
			{
				throw py::stop_iteration();
			}
			HEHandle result = _current;
			this->operator++();
			return *result;
		}

		// Begin and end methods for range-based for loop
		FaceCirculator begin() { return *this; }
		FaceCirculator end() { return FaceCirculator(); }

	private:
		HEHandle _start;
		HEHandle _current;
		bool _isEnd;
	};

	// Constant Face Circulator
	class FaceCCirculator
	{

	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = HalfEdge;
		using difference_type = std::ptrdiff_t;
		using pointer = const HalfEdge *;
		using reference = const HalfEdge &;

		FaceCCirculator() : _start{},
							_current{},
							_isEnd{true}
		{
		}
		explicit FaceCCirculator(HECHandle he) : _start{he},
														   _current{he},
														   _isEnd{false}

		{
		}
		
		FaceCCirculator &operator++()
		{
			_current = _current->next();
			if (_current == _start)
			{
				_isEnd = true; // Completed full circle
			}
			return *this;
		}

		FaceCCirculator operator++(int)
		{
			FaceCCirculator temp = *this;
			++(*this);
			// return *temp;
			return temp;
		}

		reference operator*() const
		{
			return *_current;
		}

		pointer operator->()
		{
			return &(*_current);
		}

		bool operator==(const FaceCCirculator &other) const
		{
			return (_isEnd && other._isEnd) || (_current == other._current);
		}

		bool operator!=(const FaceCCirculator &other) const
		{
			return !(*this == other);
		}

		// Begin and end methods for range-based for loop
		FaceCCirculator begin() { return *this; }
		FaceCCirculator end() { return FaceCCirculator(); }

	private:
		HECHandle _start;
		HECHandle _current;
		bool _isEnd;
	};

}


#endif