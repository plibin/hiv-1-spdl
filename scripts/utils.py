#TODO: is it OK to use this test, it assume normality, or a large enough dataset?
def ci95(x: pd.Series) -> tuple[float, float]:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) < 2:
        return (np.nan, np.nan)
    mean = x.mean()
    sem = stats.sem(x)
    lo, hi = stats.t.interval(0.95, len(x) - 1, loc=mean, scale=sem)
    return lo, hi
