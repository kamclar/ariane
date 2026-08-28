function ariane() {
    const frontend = window.ArianeFrontend;
    if (!frontend) {
        throw new Error("ARIANE frontend modules were not loaded");
    }

    return {
        ...frontend.coreState(),
        ...frontend.classificationState(),
        ...frontend.manual_reviewState(),
        ...frontend.batchState(),
        ...frontend.rulesState(),
        ...frontend.coreMethods,
        ...frontend.rulesMethods,
        ...frontend.graphsMethods,
        ...frontend.classificationMethods,
        ...frontend.formattersMethods,
        ...frontend.manual_reviewMethods,
        ...frontend.batchMethods,
    };
}

