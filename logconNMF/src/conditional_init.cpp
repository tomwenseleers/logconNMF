#include <RcppEigen.h>
#include <cmath>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__linux__)
#include <fstream>
#include <sstream>
#elif defined(__APPLE__)
#include <mach/mach.h>
#endif
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Rcpp;

// Lawson-Hanson active set. Solve the restricted problems on the design itself,
// so a pair of nearly collinear spectral directions does not square the condition number.
static VectorXd cone_solve(const MatrixXd &a, const VectorXd &b) {
    const int n = a.cols();
    VectorXd x = VectorXd::Zero(n);
    if (n == 0 || a.rows() == 0)
        return x;
    std::vector<bool> active(n, false);
    const double tol = 10 * std::numeric_limits<double>::epsilon() * std::max(a.rows(), a.cols()) *
                       std::max(1e-80, a.cwiseAbs().colwise().sum().maxCoeff()) *
                       std::max(1e-80, b.norm());
    for (int pass = 0; pass < 3 * n + 10; ++pass) {
        VectorXd grad = a.transpose() * (b - a * x);
        int add = -1;
        double best = tol;
        for (int j = 0; j < n; ++j)
            if (!active[j] && grad[j] > best) {
                best = grad[j];
                add = j;
            }
        if (add < 0)
            return x;
        active[add] = true;
        for (int it = 0; it < n + 1; ++it) {
            std::vector<int> ids;
            for (int j = 0; j < n; ++j)
                if (active[j])
                    ids.push_back(j);
            if (ids.empty())
                break;
            MatrixXd sub(a.rows(), ids.size());
            for (int k = 0; k < (int)ids.size(); ++k)
                sub.col(k) = a.col(ids[k]);
            VectorXd sol = sub.colPivHouseholderQr().solve(b);
            double step = 1;
            int drop = -1;
            for (int k = 0; k < (int)ids.size(); ++k)
                if (sol[k] <= 0) {
                    int j = ids[k];
                    double s = x[j] > 0 ? x[j] / (x[j] - sol[k]) : 0;
                    if (s <= step) {
                        step = s;
                        drop = j;
                    }
                }
            for (int k = 0; k < (int)ids.size(); ++k)
                x[ids[k]] += step * (sol[k] - x[ids[k]]);
            if (drop < 0)
                break;
            x[drop] = 0;
            active[drop] = false;
        }
    }
    return x.cwiseMax(0.);
}
// [[Rcpp::export]]
NumericVector cc_cone(const Eigen::Map<Eigen::MatrixXd> &design,
                      const Eigen::Map<Eigen::VectorXd> &y) {
    if (design.rows() != y.size())
        stop("Incompatible spectral cone.");
    return wrap(cone_solve(design, y));
}
// [[Rcpp::export]]
NumericMatrix cc_scan_unmix(const Eigen::Map<Eigen::MatrixXd> &y,
                            const Eigen::Map<Eigen::MatrixXd> &offset,
                            const Eigen::Map<Eigen::MatrixXd> &p, NumericMatrix start,
                            LogicalMatrix eligible, double phi, bool exact = false) {
    NumericMatrix out = clone(start);
    if (y.rows() != out.nrow() || p.rows() != out.ncol() || y.cols() != p.cols() ||
        offset.rows() != y.rows() || offset.cols() != y.cols() || eligible.nrow() != out.nrow() ||
        eligible.ncol() != out.ncol())
        stop("Incompatible scan regression.");
    for (int r = 0; r < out.nrow(); ++r) {
        std::vector<int> use;
        for (int k = 0; k < out.ncol(); ++k)
            if (eligible(r, k))
                use.push_back(k);
        if (use.empty())
            continue;
        const int K = use.size();
        MatrixXd a(y.cols(), K);
        VectorXd x(K);
        for (int k = 0; k < K; ++k) {
            a.col(k) = p.row(use[k]);
            x[k] = out(r, use[k]);
        }
        VectorXd yy = y.row(r), bg = offset.row(r), target = yy - bg;
        for (int it = 0; it < (exact ? 3 : 2); ++it) {
            VectorXd pred = bg + a * x;
            VectorXd w = ((!exact && it == 0) ? yy : pred)
                             .cwiseMax(exact ? 1. : std::max(1., phi))
                             .cwiseInverse();
            if (exact) {
                VectorXd root = w.cwiseSqrt();
                x = cone_solve(root.asDiagonal() * a, root.cwiseProduct(target));
                continue;
            }
            MatrixXd aw = a.transpose() * w.asDiagonal(), G = aw * a;
            VectorXd score = aw * target - G * x;
            for (int sweep = 0; sweep < 32; ++sweep) {
                double largest = 0;
                for (int k = 0; k < K; ++k) {
                    double change = std::max(0., x[k] + score[k] / std::max(G(k, k), 1e-60)) - x[k];
                    x[k] += change;
                    score.noalias() -= G.col(k) * change;
                    largest = std::max(largest, std::abs(change));
                }
                if (largest < 1e-6 * std::max(1., x.maxCoeff()))
                    break;
            }
        }
        for (int k = 0; k < K; ++k)
            out(r, use[k]) = x[k];
    }
    return out;
}
// Cyclic nonnegative least squares. Every neighbour of a coordinate occurs in
// one sorted CSC column. Ions are independent right-hand sides within each row.
// [[Rcpp::export]]
Eigen::MatrixXd cc_gram_nnls(const Eigen::MappedSparseMatrix<double> &gram,
                             const Eigen::Map<Eigen::MatrixXd> &score,
                             const Eigen::Map<Eigen::MatrixXd> &start) {
    if (gram.rows() != gram.cols() || gram.rows() != start.rows() || score.rows() != start.rows() ||
        score.cols() != start.cols())
        Rcpp::stop("Incompatible grouped-to-original spectrum regression dimensions.");
    using MatrixRows = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    MatrixRows b = start;
    if (b.size() == 0)
        return b;
    MatrixRows residual = MatrixRows::Zero(score.rows(), score.cols());
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(gram.rows());
    for (int j = 0; j < gram.outerSize(); ++j) {
        for (Eigen::MappedSparseMatrix<double>::InnerIterator a(gram, j); a; ++a) {
            if (a.row() == j)
                diagonal[j] = a.value();
            residual.row(a.row()) += a.value() * b.row(j);
        }
    }
    residual = score - residual;
    Eigen::RowVectorXd d(start.cols());
    for (int sweep = 0; sweep < 120; ++sweep) {
        double largest = 0.;
        for (int j = 0; j < b.rows(); ++j) {
            if (diagonal[j] <= 0.)
                continue;
            d = (b.row(j) + residual.row(j) / diagonal[j]).cwiseMax(0.) - b.row(j);
            const double change = d.cwiseAbs().maxCoeff();
            if (change == 0.)
                continue;
            largest = std::max(largest, change);
            b.row(j) += d;
            for (Eigen::MappedSparseMatrix<double>::InnerIterator a(gram, j); a; ++a)
                residual.row(a.row()) -= a.value() * d;
        }
        if (largest < 1e-8 * std::max(1., b.maxCoeff()))
            break;
    }
    return b;
}
// [[Rcpp::export]]
double cc_available_memory() {
#ifdef _WIN32
    MEMORYSTATUSEX m{};
    m.dwLength = sizeof(m);
    if (GlobalMemoryStatusEx(&m))
        return static_cast<double>(std::min(m.ullAvailPhys, m.ullAvailPageFile));
#elif defined(__linux__)
    std::ifstream f("/proc/meminfo");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("MemAvailable:", 0) == 0) {
            std::istringstream stream(line.substr(13));
            double kb;
            if (stream >> kb)
                return kb * 1024;
        }
#elif defined(__APPLE__)
    vm_statistics64_data_t vm{};
    mach_msg_type_number_t n = HOST_VM_INFO64_COUNT;
    mach_port_t host = mach_host_self();
    vm_size_t page = 0;
    if (host_page_size(host, &page) == KERN_SUCCESS &&
        host_statistics64(host, HOST_VM_INFO64, reinterpret_cast<host_info64_t>(&vm), &n) ==
            KERN_SUCCESS)
        return static_cast<double>(page) * (vm.free_count + vm.inactive_count);
#endif
    return NA_REAL;
}
