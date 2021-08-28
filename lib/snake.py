from collections import defaultdict
from warnings import warn
from numpy import ceil

alias_recipe = "ln -rs {input} {output}"
alias_recipe_norelative = "ln -sT {input} {output}"
alias_fmt = lambda input, output: alias_recipe.format(input=input, output=output)
curl_recipe = "curl '{params.url}' > {output}"
curl_unzip_recipe = "curl '{params.url}' | zcat > {output}"

limit_numpy_procs = """
        export MKL_NUM_THREADS={threads}
        export OPENBLAS_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        export VECLIB_MAXIMUM_THREADS={threads}

        """
limit_numpy_procs_to_1 = limit_numpy_procs.format(threads=1)

independent_theano_compiledir = """
        # Circumvent theano compiledir locking.
        compiledir=$(mktemp -d)
        # tar --strip-components 1 -xzf raw/theano.tgz -C $compiledir
        export THEANO_FLAGS="base_compiledir=$compiledir"
        echo $THEANO_FLAGS

        """

# Utility wildcard constrains
noperiod_wc = "[^.]+"
integer_wc = "[0-9]+"
float_noperiod_wc = "[0-9]+(e[0-9]+)?"
single_param_wc = "[^.-]+"
params_wc = noperiod_wc


def resource_calculator(
    baseline=1,
    threads_exponent=0,
    attempt_base=1,
    input_size_exponent=None,
    **input_size_multipliers
):
    if input_size_exponent is None:
        input_size_exponent = {}
    _input_size_exponent = defaultdict(lambda: 1)
    _input_size_exponent.update(input_size_exponent)
    input_size_exponent = _input_size_exponent

    def func(wildcards, input, threads, attempt):
        input_sizes = {}
        for k in input_size_multipliers:
            input_sizes[k] = getattr(input, k).size / 1024 / 1024
            # # print(type(getattr(input, k).size))
            # try:
            # except AttributeError as err:
            #     warn(str(err))
            #     input_sizes[k] = 0
        base_estimate = baseline + sum(
            [input_sizes[k] * input_size_multipliers[k] for k in input_size_multipliers]
        )
        outvalue = base_estimate * threads ** threads_exponent * attempt_base ** attempt
        # print(weighted_input_size, threads, threads_exponent, attempt_base, attempt, outvalue)
        return ceil(outvalue)

    return func
