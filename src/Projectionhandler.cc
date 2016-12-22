#include "lsst/jointcal/Projectionhandler.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/CcdImage.h"

namespace lsst {
namespace jointcal {


class Mapping;


/**********   Stuff for providing Sk22TP gtransfos to a DistortionModel ***/

OneTPPerVisitHandler::OneTPPerVisitHandler(const CcdImageList &ccdImageList)
{
  for (auto i=ccdImageList.cbegin(); i!= ccdImageList.end(); ++i)
    {
      const CcdImage &im = **i;
      if (tMap.find(im.getVisit()) == tMap.end())
	tMap[im.getVisit()] = im.Sky2TP()->Clone();
    }

}

const Gtransfo* OneTPPerVisitHandler::Sky2TP(const CcdImage &ccdImage) const
{
  auto it=tMap.find(ccdImage.getVisit());
  if (it==tMap.end()) return nullptr;
  return &*(it->second);
}

}} // end of namespaces
