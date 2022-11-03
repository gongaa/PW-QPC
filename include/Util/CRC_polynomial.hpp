#ifndef CRC_POLYNOMIAL_HPP_
#define CRC_POLYNOMIAL_HPP_

#include <string>
#include <vector>
#include <tuple>
#include <map>

template <typename B = int>
class CRC_polynomial
{
protected:
	const static std::map<std::string, std::tuple<unsigned, int>> known_polynomials;
	std::vector<B> polynomial;
	unsigned       polynomial_packed;
	std::vector<B> buff_crc;
	const int size;
	const int K;

public:
	CRC_polynomial(const int K, const std::string &poly_key, const int size = 0);
	virtual ~CRC_polynomial() = default;

	static int         get_size (const std::string &poly_key);
	static std::string get_name (const std::string &poly_key);
	static unsigned    get_value(const std::string &poly_key);

	virtual void build       (const B *U_K1, B *U_K2, const size_t frame_id);
	virtual void extract     (const B *V_K1, B *V_K2, const size_t frame_id);
	virtual bool check       (const B *V_K          , const size_t frame_id);
	// virtual bool check_packed(const B *V_K          , const size_t frame_id);

	void generate(const B *U_in, B *U_out,
	              const int off_in, const int off_out,
	              const int loop_size);
};

#include "CRC_polynomial.hxx"

#endif /* CRC_POLYNOMIAL_HPP_ */
