// Generated by rstantools.  Do not edit by hand.

/*
    test is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    test is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with test.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_logn_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_logn");
    reader.add_event(37, 35, "end", "model_logn");
    return reader;
}
#include <stan_meta_header.hpp>
class model_logn
  : public stan::model::model_base_crtp<model_logn> {
private:
        int N;
        vector_d y;
        double L_m;
        double U_m;
        double L_cv;
        double U_cv;
public:
    model_logn(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_logn(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_logn_namespace::model_logn";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 3;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(N));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "L_m", "double", context__.to_vec());
            L_m = double(0);
            vals_r__ = context__.vals_r("L_m");
            pos__ = 0;
            L_m = vals_r__[pos__++];
            check_greater_or_equal(function__, "L_m", L_m, 0);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "U_m", "double", context__.to_vec());
            U_m = double(0);
            vals_r__ = context__.vals_r("U_m");
            pos__ = 0;
            U_m = vals_r__[pos__++];
            check_greater_or_equal(function__, "U_m", U_m, 0);
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "L_cv", "double", context__.to_vec());
            L_cv = double(0);
            vals_r__ = context__.vals_r("L_cv");
            pos__ = 0;
            L_cv = vals_r__[pos__++];
            check_greater_or_equal(function__, "L_cv", L_cv, 0);
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "U_cv", "double", context__.to_vec());
            U_cv = double(0);
            vals_r__ = context__.vals_r("U_cv");
            pos__ = 0;
            U_cv = vals_r__[pos__++];
            check_greater_or_equal(function__, "U_cv", U_cv, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 11;
            num_params_r__ += 1;
            current_statement_begin__ = 12;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_logn() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 11;
        if (!(context__.contains_r("m")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable m missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("m");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "m", "double", context__.to_vec());
        double m(0);
        m = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, m);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable m: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 12;
        if (!(context__.contains_r("cv")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable cv missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("cv");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "cv", "double", context__.to_vec());
        double cv(0);
        cv = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, cv);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable cv: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 11;
            local_scalar_t__ m;
            (void) m;  // dummy to suppress unused var warning
            if (jacobian__)
                m = in__.scalar_lb_constrain(0, lp__);
            else
                m = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 12;
            local_scalar_t__ cv;
            (void) cv;  // dummy to suppress unused var warning
            if (jacobian__)
                cv = in__.scalar_lb_constrain(0, lp__);
            else
                cv = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 16;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 17;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            stan::math::initialize(sigma, DUMMY_VAR__);
            stan::math::fill(sigma, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 19;
            stan::math::assign(mu, stan::math::log((m / stan::math::sqrt((1 + pow(cv, 2))))));
            current_statement_begin__ = 20;
            stan::math::assign(sigma, stan::math::sqrt(stan::math::log((1 + pow(cv, 2)))));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 16;
            if (stan::math::is_uninitialized(mu)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: mu";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            current_statement_begin__ = 17;
            if (stan::math::is_uninitialized(sigma)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: sigma";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable sigma: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            check_greater_or_equal(function__, "sigma", sigma, 0);
            // model body
            current_statement_begin__ = 24;
            lp_accum__.add(uniform_log<propto__>(m, L_m, U_m));
            current_statement_begin__ = 25;
            lp_accum__.add(uniform_log<propto__>(cv, L_cv, U_cv));
            current_statement_begin__ = 26;
            lp_accum__.add(normal_log<propto__>(stan::math::log(y), mu, sigma));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("m");
        names__.push_back("cv");
        names__.push_back("mu");
        names__.push_back("sigma");
        names__.push_back("x_tilde");
        names__.push_back("y_tilde");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_logn_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double m = in__.scalar_lb_constrain(0);
        vars__.push_back(m);
        double cv = in__.scalar_lb_constrain(0);
        vars__.push_back(cv);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 16;
            double mu;
            (void) mu;  // dummy to suppress unused var warning
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 17;
            double sigma;
            (void) sigma;  // dummy to suppress unused var warning
            stan::math::initialize(sigma, DUMMY_VAR__);
            stan::math::fill(sigma, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 19;
            stan::math::assign(mu, stan::math::log((m / stan::math::sqrt((1 + pow(cv, 2))))));
            current_statement_begin__ = 20;
            stan::math::assign(sigma, stan::math::sqrt(stan::math::log((1 + pow(cv, 2)))));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 17;
            check_greater_or_equal(function__, "sigma", sigma, 0);
            // write transformed parameters
            if (include_tparams__) {
                vars__.push_back(mu);
                vars__.push_back(sigma);
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 30;
            double x_tilde;
            (void) x_tilde;  // dummy to suppress unused var warning
            stan::math::initialize(x_tilde, DUMMY_VAR__);
            stan::math::fill(x_tilde, DUMMY_VAR__);
            current_statement_begin__ = 31;
            double y_tilde;
            (void) y_tilde;  // dummy to suppress unused var warning
            stan::math::initialize(y_tilde, DUMMY_VAR__);
            stan::math::fill(y_tilde, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 33;
            stan::math::assign(x_tilde, normal_rng(mu, sigma, base_rng__));
            current_statement_begin__ = 34;
            stan::math::assign(y_tilde, stan::math::exp(x_tilde));
            // validate, write generated quantities
            current_statement_begin__ = 30;
            vars__.push_back(x_tilde);
            current_statement_begin__ = 31;
            check_greater_or_equal(function__, "y_tilde", y_tilde, 0);
            vars__.push_back(y_tilde);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_logn";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "m";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "cv";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "x_tilde";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "y_tilde";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "m";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "cv";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "x_tilde";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "y_tilde";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_logn_namespace::model_logn stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
