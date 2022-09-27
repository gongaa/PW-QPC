#ifndef ENCODER_RM_CSS_HPP_
#define ENCODER_RM_CSS_HPP_

class Encoder_RM_CSS 
{
protected:
    int m, rx, rz;
    int N, K;

public:
    Encoder_RM_CSS(const int& m, const int& rx, const int& rz);
    int get_K() { return K; }
    int get_N() { return N; }
    bool is_logical_X(const int *X_N);
    bool is_X_stabilizer(const int *X_N);
    bool is_logical_Z(const int *X_N);
    bool is_Z_stabilizer(const int *X_N);
};

#endif // ENCODER_RM_CSS_HPP_