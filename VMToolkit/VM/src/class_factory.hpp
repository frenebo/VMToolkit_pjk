/*!
 * \file class_factory.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief A simple class factory
 */


#ifndef __FACTORY_TYPES_HPP__
#define __FACTORY_TYPES_HPP__

#include <memory>
#include <map>
#include <string>
#include <stdexcept>

/*

	A simple implementation of a class factory.
	The assumption is that one has a class hierarchy with a
	base class (BaseType) and its children (DerivedType).
	Class factory then creates a mapping (map) between a string
	(class name) and a smart pointer to an object of a give base or
	derived class.

	Usage:

	   typedef ClassFactory<Base> class_factory;

	   void f()
	   {
		 // ...
		 class_factory cf;

		 cf.add<Base, int>("base", 1);
		 cf.add<ChildA, int, int>("child", 1, 2);

		 // ...

		 for (auto obj : cf.factory())
			obj.second->...;   // ... stands for the member function call
	   }

*/
namespace VMTutorial
{
	template <typename BaseType>
	class ClassFactory
	{

		typedef std::unique_ptr<BaseType> ptrBaseType;
		typedef std::map<std::string, ptrBaseType> factory_type;

	public:
		ClassFactory() = default;
		~ClassFactory() = default;

		// Adds a new object to the factory. The constructor of the DerivedType is invoked.
		template <typename DerivedType, typename... Args>
		void add(const std::string &key, const Args &...args);

		// Remove
		void remove(const std::string &key)
		{
			if (factory_map.find(key) != factory_map.end())
				factory_map.erase(key);
			else
				throw std::invalid_argument(key + " is not in the class factory.");
		}

		// Retrun pointer to the object associated with the key
		ptrBaseType &get(const std::string &key)
		{
			if (factory_map.find(key) != factory_map.end())
				return factory_map[key];
			else
				throw std::invalid_argument(key + " is not in the class factory.");
		}

		// Return the entire factory
		const factory_type &factory() const { return factory_map; }

	protected:
		factory_type factory_map;
	};

	template <typename BaseType>
	template <typename DerivedType, typename... Args>
	void ClassFactory<BaseType>::add(const std::string &key, const Args &...args)
	{
		if (factory_map.find(key) != factory_map.end()) // Add only if the key is not already in the class factory.
			throw std::invalid_argument(key + " is already in the class factory.");
		factory_map[key] = std::make_unique<DerivedType>(args...);
	}
}
#endif
