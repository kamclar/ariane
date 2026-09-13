function ariane() {
    const frontend = window.ArianeFrontend;
    if (!frontend) {
        throw new Error("ARIANE frontend modules were not loaded");
    }

    return frontend.composeStrict([
        { name: "core state", value: frontend.coreState() },
        { name: "classification state", value: frontend.classificationState() },
        { name: "manual-review state", value: frontend.manual_reviewState() },
        { name: "batch state", value: frontend.batchState() },
        { name: "rules state", value: frontend.rulesState() },
        { name: "core methods", value: frontend.coreMethods },
        { name: "rules methods", value: frontend.rulesMethods },
        { name: "graph methods", value: frontend.graphsMethods },
        { name: "classification methods", value: frontend.classificationMethods },
        { name: "formatter methods", value: frontend.formattersMethods },
        { name: "manual-review methods", value: frontend.manual_reviewMethods },
        { name: "manual-review API methods", value: frontend.manualReviewApiMethods },
        { name: "manual-review PS1 methods", value: frontend.manualReviewPs1Methods },
        { name: "manual-review persistence methods", value: frontend.manualReviewPersistenceMethods },
        { name: "batch methods", value: frontend.batchMethods },
    ]);
}
