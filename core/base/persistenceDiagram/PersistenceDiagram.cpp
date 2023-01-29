#include <PersistenceDiagram.h>

using namespace std;
using namespace ttk;

PersistenceDiagram::PersistenceDiagram() {
  setDebugMsgPrefix("PersistenceDiagram");
}

void ttk::PersistenceDiagram::sortPersistenceDiagram(
  std::vector<PersistencePair> &diagram, const SimplexId *const offsets) const {

  auto cmp = [offsets](const PersistencePair &a, const PersistencePair &b) {
    return offsets[a.birth.id] < offsets[b.birth.id];
  };

  std::sort(diagram.begin(), diagram.end(), cmp);
}
