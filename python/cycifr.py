def get_gmm_and_pos_label(
    array, n_components=2, n_steps=5000
):
    """
    Compute the Gaussian Mixture Model (GMM) fitting of a given dataset and
    identify the positive class based on the GMM means.

    Parameters
    ----------
    array : numpy array
        The input dataset to be modeled as a GMM.
    n_components : int, optional (default=2)
        The number of components to use in the GMM.
    n_steps : int, optional (default=5000)
        The number of points to use to sample the GMM PDF.

    Returns
    -------
    gmm : GaussianMixture object
        A trained GMM object from the scikit-learn library.
    label : int
        The positive class label identified based on the GMM means.
    cutoffs : numpy array
        An array of cutoff points that can be used to threshold the data.

    Notes
    -----
    This function fits a GMM to the input data using the `sklearn.mixture.GaussianMixture`
    object. It then identifies the positive class based on the GMM means and computes
    a set of cutoff points that can be used to threshold the data.

    The `n_components` and `n_steps` parameters can be used to control the complexity of
    the GMM and the number of points used to sample the PDF.

    """
    # Fit a GMM to the input data
    gmm = sklearn.mixture.GaussianMixture(n_components=2, covariance_type='spherical', random_state=0)
    gmm.fit(array.reshape(-1, 1))
    gmm$weights_

    # Identify the positive class based on the GMM means
    label = np.argmax(gmm.means_)

    # Compute a set of cutoff points that can be used to threshold the data
    low = gmm.means_.min() - 2*np.sqrt(gmm.covariances_[np.argmin(gmm.means_)])
    high = gmm.means_.max() + 2*np.sqrt(gmm.covariances_[np.argmax(gmm.means_)])
    ref_space = np.linspace(low, high, n_steps)
    result = gmm.predict(ref_space.reshape(-1, 1))
    idx = np.where(np.ediff1d(result) != 0)
    cutoffs = ref_space[idx]

    # Return the GMM, positive label, and cutoffs
    return gmm, label, cutoffs


def plot_gmm_fitting(array, gmm, ax):
    """
    Plot the Gaussian Mixture Model (GMM) fitting of a given dataset.

    Parameters
    ----------
    array : numpy array
        The input dataset to be modeled as a GMM.
    gmm : GaussianMixture object
        A trained GMM object from the scikit-learn library.
    ax : matplotlib Axes object
        The target axes to plot the GMM fitting.

    Returns
    -------
    ax : matplotlib Axes object
        The same input axes object with the GMM fitting plotted.

    Notes
    -----
    This function plots the histogram of the input dataset, and the GMM fitting
    for the data. The GMM fitting is computed by evaluating the probability
    density function (PDF) of the GMM at a set of equally spaced points, and
    then plotting the individual PDFs of each component in the GMM, as well
    as the overall PDF of the GMM.

    """
    # Set the current axes to the target axes
    plt.sca(ax)

    # Plot the histogram of the input data
    _ = plt.hist(array.flatten(), color='lightgray', bins=200, density=True)

    # Evaluate the PDF of the GMM at a set of points
    x = np.linspace(array.min(), array.max(), 200)
    log_prob = gmm.score_samples(x.reshape(-1, 1))
    responsibilities = gmm.predict_proba(x.reshape(-1, 1))
    pdf = np.exp(log_prob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]

    # Plot the individual PDFs of each component in the GMM
    mean_index = np.argmax(pdf_individual, axis=0)
    rank_map = mean_index.argsort().argsort()
    ax.set_prop_cycle(color=plt.get_cmap('Dark2')(rank_map))
    ax.plot(x, pdf_individual)

    # Plot the overall PDF of the GMM
    ax.plot(x, pdf, '--k')

    # Return the target axes with the GMM fitting plotted
    return ax
  
def plot_hist_gmm(
    df,
    markers,
    n_components=2,
    subplot_grid_shape=None,
    transform_log=True,
    xlim_percentiles=(0, 100),
    cum_density=False,
    hide_yaxis_left=True
):  
    if transform_log:
        df = df.transform(np.log1p)
        revert_func = np.expm1
    else:
        revert_func = np.array
    if subplot_grid_shape is None:
        subplot_grid_shape = (1, len(markers))
    n_rows, n_cols = subplot_grid_shape
    fig, axes = plt.subplots(n_rows, n_cols, sharex=True)
    axes = np.array(axes)

    for m, ax in zip(markers, axes.ravel()):
        gmm, _, cutoffs = get_gmm_and_pos_label(
            df[m].values, n_components=n_components
        )
        plot_gmm_fitting(df[m].values, gmm, ax)
        ax.title.set_text(m)
        if hide_yaxis_left:
            ax.yaxis.set_visible(False)

        p1, p2 = np.array(xlim_percentiles) / 100
        axis_min = df.loc[:, markers].quantile(p1).min()
        axis_max = df.loc[:, markers].quantile(p2).max()

        color_cum = 'gray'

        pax = ax.twinx()
        pax = plot_cumulative(
            df[m].values, pax, 
            hist_kwargs=dict(color=color_cum, density=cum_density)
        )
        pax.tick_params(axis='y', labelsize=8, colors=color_cum)

        print(cutoffs)

        cutoff_range = np.ptp(cutoffs)
        if cutoff_range == 0: cutoff_range = 1
        cutoff_colors = plt.get_cmap('plasma')(
            (cutoffs - np.min(cutoffs)) / cutoff_range
        )

        for co, cc in zip(cutoffs, cutoff_colors):
            ax.axvline(x=co, c=cc, alpha=0.2)
            ax.annotate(
                '',
                xy=(co, 0), xytext=(co, -0.05),
                xycoords=('data', 'axes fraction'),
                arrowprops=dict(arrowstyle='wedge, tail_width=0.7, shrink_factor=0.5', color=cc)
            )
        ax.set_xlim(axis_min, axis_max)
        # cutoff_string = np.round(revert_func(cutoffs)).astype(int)

        for i, (co, cc) in enumerate(
            zip(revert_func(cutoffs)[::-1], cutoff_colors[::-1])
        ):
            text = ax.text(
                ax.get_xlim()[0] + 0.02*np.diff(ax.get_xlim()), 
                ax.get_ylim()[1] - 0.05*(i+1)*np.diff(ax.get_ylim()), 
                f'{np.round(co).astype(int)}', 
                fontsize=10, c=cc
            )
            text_outline = mpatheffects.Stroke(linewidth=1, foreground='#000')
            text.set_path_effects(
                [text_outline, mpatheffects.Normal()]
            )
    plt.tight_layout()
    for aax in fig.axes:
        aax.spines['right'].set_color(color_cum)
        power_label = aax.yaxis.get_offset_text()
        power_label.set_visible(False)
        aax.annotate(
            power_label.get_text(), xy=(1.02, 1.01),
            xycoords='axes fraction', fontsize=10,
            color=color_cum
        )
    plt.sca(ax)
