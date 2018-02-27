#pragma once

/*********************************************
Framework to test all functionalities,
defined in BinaryTrees.cpp and functions.cpp,
for working with compressed matrices
(for example, Add, Mult, Inverse and etc.)
*********************************************/

#include "definitions.h"

using namespace std;

// Interface

template <class T>
ostream& operator << (ostream& os, const vector<T>& s);

template <class T>
ostream& operator << (ostream& os, const set<T>& s);

template <class K, class V>
ostream& operator << (ostream& os, const map<K, V>& m);

template<class T, class U>
void AssertEqual(const T& t, const U& u, const string& hint = {});

template<class T, class U>
void AssertLess(const T& t, const U& u, const string& hint = {});

void Assert(bool b, const string& hint);


class TestRunner {
public:

	template <class TestFunc, class T>
	void RunOld(TestFunc(*func)(T&), const string& test_name);

	template <class Shell, class One>
	void RunTest(Shell shell, One oneTest, const string& test_name);

	template <class Shell, class Simple, class Struct>
	void RunTests(Shell shell, Simple testSimple, Struct testStruct, const string& test_name);

	int GetAll() const;
	int GetPassed() const;
	int GetFailed() const;

	~TestRunner();

private:
	int all_tests = 0;
	int passed = 0;
	int failed = 0;
};


// Implementation

template <class T>
ostream& operator << (ostream& os, const vector<T>& s) {
	os << "{";
	bool first = true;
	for (const auto& x : s) {
		if (!first) {
			os << ", ";
		}
		first = false;
		os << x;
	}
	return os << "}";
}

template <class T>
ostream& operator << (ostream& os, const set<T>& s) {
	os << "{";
	bool first = true;
	for (const auto& x : s) {
		if (!first) {
			os << ", ";
		}
		first = false;
		os << x;
	}
	return os << "}";
}

template <class K, class V>
ostream& operator << (ostream& os, const map<K, V>& m) {
	os << "{";
	bool first = true;
	for (const auto& kv : m) {
		if (!first) {
			os << ", ";
		}
		first = false;
		os << kv.first << ": " << kv.second;
	}
	return os << "}";
}

template<class T, class U>
void AssertEqual(const T& t, const U& u, const string& hint) {
	if (t != u) {
		ostringstream os;
		os << "Assertion failed: " << t << " != " << u;
		if (!hint.empty()) {
			os << " hint: " << hint;
		}
		throw runtime_error(os.str());
	}
}

template<class T, class U>
void AssertLess(const T& t, const U& u, const string& hint) {
	if (t > u) {
		ostringstream os;
		os << "FAILED. "  << t << " > " << u;
		if (!hint.empty()) {
			os << " hint: " << hint;
		}
		throw runtime_error(os.str());
	}
}

void Assert(bool b, const string& hint) {
	AssertEqual(b, true, hint);
}


template <class TestFunc, class T>
void TestRunner::RunOld(TestFunc(*func)(T&), const string& test_name) {
	
	T n = 0;
	try 
	{
		(*func)(n);
	}
	catch (runtime_error& e)
	{
		++fail_count;
		cerr << test_name << " fail: " << e.what() << endl;
	}
	catch (...) {
		++fail_count;
		cerr << "Unknown exception caught" << endl;
	}
	cerr << test_name << ". " << n << " tests passed. OK" << endl;
}

template <class Shell, class One>
void TestRunner::RunTest(Shell shell, One oneTest, const string& test_name)
{
	int local_tests = 0;
	int local_fails = 0;
	shell(oneTest, test_name, local_tests, local_fails);

	all_tests += local_tests;
	failed += local_fails;
	passed = all_tests - failed;

	cerr << test_name << ". Tests: " << local_tests << ". Passed: " << local_tests - local_fails;
	cerr << ". Failed: " << local_fails << ".";
	if (local_fails == 0) cerr << " OK.";
	cerr << endl;

}

template <class Shell, class Simple, class Struct>
void TestRunner::RunTests(Shell shell, Simple testSimple, Struct testStruct, const string& test_name)
{
	int local_tests = 0;
	int local_fails = 0;
	shell(testSimple, test_name, local_tests, local_fails);
	shell(testStruct, test_name, local_tests, local_fails);

	all_tests += local_tests;
	failed += local_fails;
	passed = all_tests - failed;

	cerr << test_name << ". Tests: " << local_tests << ". Passed: " << local_tests - local_fails;
	cerr << ". Failed: " << local_fails << ".";
	if (local_fails == 0) cerr << " OK.";
	cerr << endl;
}

TestRunner::~TestRunner()
{
	if (failed > 0)
	{
		cerr << failed << " unit tests failed. Terminate" << endl;
		system("pause");
		exit(1); // окончание работы программы, если количество упавших тестов > 0
	}
}

int TestRunner::GetAll() const
{
	return all_tests;
}

int TestRunner::GetPassed() const
{
	return passed;
}

int TestRunner::GetFailed() const
{
	return failed;
}