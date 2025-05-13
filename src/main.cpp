#include "matrix.h"
#include "vector.h"
using namespace supp_math;
template<typename T, typename U>
using ParamType = std::tuple<Matrix<T>, Matrix<U>, Matrix<std::common_type_t<T, U>>>;

template<typename T, typename U>
ParamType<T, U> get_sample_tuple()
{

	return std::make_tuple(
	    Matrix<T>({
	        {0, 2, 3},
          {3, 8, 7},
          {6, 6, 6}
  }),
	    Matrix<U>({{5, 6}, {7, 8}}), Matrix<std::common_type_t<T, U>>({{6, 8}, {10, 12}})
	);
}

int main()
{
	// auto a = get_sample_tuple<int, int>();

	// Matrix<double> mmmmmmm = std::get<0>(a);
	supmath::Vector<double> v1     = supmath::Vector<double>({1, 2, 3});
	supmath::Vector<int>    v2     = supmath::Vector<int>({1, 5, 7});
	supmath::Vector         v3     = supmath::Vector(v2, 1);

	auto                    rrrrrr = v1.projectionOnto(v2);
	std::cout << v3 << std::endl;

	supp_math::Matrix<double> sample_double_4x4 = supp_math::Matrix<double>({
	    {  6,   -4, 2.11, 6},
	    { -7, 6.88,  7.1, 1},
	    {  6,   -4,    8, 3},
	    {9.6,  -11,    8, 9},
	});

	auto                      rrrdsklvsd        = sample_double_4x4 * sample_double_4x4;
	// printf("%f", sample_double_4x4.determinant());
	// Matrix<double> sample_double_2x3 = Matrix<double>({
	//           {  6.3, 8, -4.6},
	//           {-7.77, 9,    1}
	//       });
	//       sample_double_2x3.determinant();
	//   std::cout << mmmmmmm << std::endl;

	// auto inv = mmmmmmm.inverse();
	// std::cout << inv << std::endl;
	// std::cout << mmmmmmm.determinant() << std::endl;

	// auto [U, L, P] = gaussian_elimination(mmmmmmm);
	// std::cout << U << std::endl;
	// std::cout<< L << std::endl;
	// std::cout << P << std::endl;
	// Matrix<double> sample_double_7x7 = Matrix<double>({
	//     {    6,    -4,  2.11,     6,   9.3,     4,     5},
	//     { 8.88,   9.3,  10.1,   2.5,  -6.3, -15.4, -0.23},
	//     {    6,    -4,     8,     3,    -9,   5.5,  0.13},
	//     {  9.6,   -11,     8,     9, -0.99,  0.99,     2},
	//     { -6.2,     1,  3.33,  8.12,     6,  7.89,     3},
	//     {-8.21, -9.63,  7.77, -1.23,     6,   9.6,    -4},
	//     { 7.13,  6.54, -8.96, -4.32,  0.15,  9.88,  -4.2},
	// });
	// auto n = sample_double_7x7.inverse();
	// std::cout << n << std::endl;

	// Matrix<double> m1 = Matrix({
	//	{0, 2, 3, 4},
	//       {  0,  2,  3,  4},
	//	{6.3, 7, 8,6},
	//	{9, 10, 11, 12}
	//});

	//	std::cout<< m1.inverse() << std::endl;
	// m1.swap(0, 2);
	//       m1.transpose();

	// std::cout << m1.determinant() << std::endl;
	// auto abc = Matrix<int>({{1, 2}, {3, 4}});
	// Matrix<int> mmmm1;
	// mmmm1 = abc;

	// std::cout << &abc[0] << "   " << &mmmm1[0] << std::endl;

	// mmmm1.p();

	// std::vector<std::vector<std::int16_t>> v1 =
	//	{
	//		{1,2},
	//		{3,4}
	//	};
	// Matrix<std::int16_t> m1(v1);

	// std::vector<std::vector<std::int32_t>> v2 =
	//	{
	//		{10,20},
	//		{30,40}
	//	};
	// Matrix<std::int32_t> m2(v2);

	// auto rrr = m2 + m1;
	// rrr.p();

	// std::vector<std::vector<double>> v3 =
	//{
	//	{10.456,20.123},
	//	{30.49,40}
	// };
	// Matrix m3(v3);
	//       Matrix<double> m4;
	//       m4 = Matrix<double>({
	//           {-12, 26},
	//           {-30, 47}
	//       });
	//   m4.p();

	// auto rrrr = m3 + m1;
	// auto rrrr5 = m3 - m1;
	// auto rrrr6 = m3 * m1;
	// auto rrrr7 = m1 - m3;
	// auto rrrr8 = m1 - m3;

	// auto rrrr9 = m1 - m3;

	// auto rrrr1 = 3 * m1;
	// auto rrrr2 = 6.3 * m1;
	// auto rrrr3 = m1 * 5;
	// auto rrrr4 = m1 * 3.33;

	// Matrix rrrr10 = m1.element_wise_product(m2);
	//
	// rrrr10.p();

	// Matrix<double> m10 = m2;
	// m10[0][0] = 10.2222;
	// m10.p();
	////m3.p();
	////std::cout << "\n";
	////Matrix<int> m4 = m3;
	////m4.p();

	// std::vector<std::vector<int>> v5 =
	//{
	//	{3,4,6},
	//	{1,2,7}
	// };
	// Matrix<int> m5(v5);
	// auto aassada = m1 + m5;
	// std::vector<std::vector<float>> v6 =
	//{
	//	{3.9,4,6.5},
	//	{1.4,2.5,7},
	//	{4,3,1.4},
	// };
	// Matrix<float> m6(v6);

	// auto rr = m5.transpose();
	// rr.p();
	// auto re = m5 * m6;
	// re.p();

	// float x = 3.2;
	// auto re2 = x * m1;
	// re2.p();

	// std::cout << re[1][0] << "\n";

	// re[1][0] = -100;

	// std::cout << re[1][0] << "\n";

	// if (&m1 == &re) {
	//	std::cout << "m1 and re are the same object" << std::endl;
	// }
	// else {
	//	std::cout << "m1 and re are not the same object" << std::endl;
	// }

	return 0;
}
