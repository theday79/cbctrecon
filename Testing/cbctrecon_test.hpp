#ifndef CBCTRECON_TEST_HPP
#define CBCTRECON_TEST_HPP

/* Defines necessary for using TinyRefl */

#define TINYREFL_API_CODEGEN_VERSION_MAJOR 0
#define TINYREFL_API_CODEGEN_VERSION_MINOR 1
#define TINYREFL_API_CODEGEN_VERSION_FIX   1

// Test backend for tinyrefl metadata
#define TINYREFL_GODMODE(...) // No Gods here
#define TINYREFL_SEQUENCE(...) __VA_ARGS__
#define TINYREFL_STRING(...) __VA_ARGS__
#define TINYREFL_TYPE(name, fullname) fullname
#define TINYREFL_VALUE(type, value) // We don't care about values (pointers to members, enums, etc)
#define TINYREFL_ATTRIBUTE(name, namespace_, full_attribute, args) // We don't care about attributes
#define TINYREFL_CONSTRUCTOR(name, fullname, class_type, signature, attributes) // We don't care about ctors
#define TINYREFL_MEMBER_FUNCTION(name, fullname, parent_class_type, return_type, signature, pointer, attributes) // we don't care about member functions
#define TINYREFL_MEMBER_VARIABLE(name, fullname, parent_class_type, value_type, pointer, attributes) (value_type, name)
#define TINYREFL_ENUM_VALUE(name, fullname, type, value, attributes) // we don't care about enums
#define TINYREFL_REFLECT_MEMBER(member) // We don't care about metadata definitions of members
#define TINYREFL_REFLECT_ENUM_VALUE(value) // we don't care about enums
#define TINYREFL_REFLECT_ENUM(name, type, values, attributes) // we don't care about enums
#define TINYREFL_REFLECT_CLASS(classname, classtype, bases, constructors, member_functions, member_variables, classes, enums, attributes) // \
//    BOOST_FUSION_ADAPT_STRUCT(classname, TINYREFL_PP_UNWRAP member_variables);

#endif // CBCTRECON_TEST_HPP
