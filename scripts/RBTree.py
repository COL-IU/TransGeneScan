#! /usr/bin/python

RED = True
BLACK = False

class RBNode:
 def __init__(self,interval,obj):
  self.color = None
  self.left = None
  self.right = None
  self.p = None
  self.key = {}
  self.key['low'] = interval[0]
  self.key['high'] = interval[1]
  self.obj = obj
  self.max = 0
 
class RBTree:
 def __init__(self,items):
  self.nil = RBNode([0,0],None)
  self.root = self.nil
  if len(items) > 0 :
   for item in items:
    tempNode = RBNode(item[0:2],item[2])
    self.insert(tempNode)
 def fixMax(self,x):
  temp = x
  while temp != None:
   leftMax = 0
   rightMax = 0
   if temp.left != None:
    leftMax = temp.left.max
   if temp.right != None:
    rightMax = temp.right.max
   temp.max = max(temp.key['high'],leftMax,rightMax) 
   temp = temp.p
 def leftRotate(self,x):
  y = x.right
  x.right = y.left
  if y.left != self.nil:
   y.left.p = x
  y.p = x.p
  if x.p == self.nil:
   self.root = y
  elif x == x.p.left:
   x.p.left = y
  else:
   x.p.right = y
  y.left = x
  x.p = y
  self.fixMax(x)
 def rightRotate(self,y):
  x = y.left
  y.left = x.right
  if x.right != self.nil:
   x.right.p = y
  x.p = y.p
  if y.p == self.nil:
   self.root = x
  elif y == y.p.left:
   y.p.left = x
  else:
   y.p.right = x
  x.right = y
  y.p = x
  self.fixMax(y)
 def insertFixup(self,z):
  while z.p.color == RED:
   if z.p == z.p.p.left:
    y = z.p.p.right
    if y.color == RED:
     z.p.color = BLACK
     y.color = BLACK
     z.p.p.color = RED
     z = z.p.p
    else:
     if z == z.p.right:
      z = z.p
      self.leftRotate(z)
     z.p.color = BLACK
     z.p.p.color = RED
     self.rightRotate(z.p.p)
   else:
    y = z.p.p.left
    if y.color == RED:
     z.p.color = BLACK
     y.color = BLACK
     z.p.p.color = RED
     z = z.p.p
    else:
     if z == z.p.left:
      z = z.p
      self.rightRotate(z)
     z.p.color = BLACK
     z.p.p.color = RED
     self.leftRotate(z.p.p)
  self.root.color = BLACK

 def insert(self,z):
  y = self.nil
  x = self.root
  while x != self.nil:
   y = x
   if z.key['low'] < x.key['low']:
    x = x.left
   else:
    x = x.right
  z.p = y
  if y == self.nil:
   self.root = z
  elif z.key < y.key:
   y.left = z
  else:
   y.right = z
  z.left = self.nil
  z.right = self.nil
  z.color = RED
  z.max = z.key['high']
  self.fixMax(z.p)
  self.insertFixup(z)

 def inOrderPrint(self,z):
  if z == self.nil:
   return
  self.inOrderPrint(z.left)
  print "(["+str(z.key['low'])+","+str(z.key['high'])+"],"+str(z.color)+","+str(z.max)+","+str(z.frame)+")"
  print "\tLeft: (["+str(z.left.key['low'])+","+str(z.left.key['high'])+"],"+str(z.left.color)+","+str(z.left.max)+","+str(z.left.frame)+")"
  print "\tRight: (["+str(z.right.key['low'])+","+str(z.right.key['high'])+"],"+str(z.right.color)+","+str(z.right.max)+","+str(z.right.frame)+")"
  self.inOrderPrint(z.right)
 def overlap(self,x,y):
  if x["low"] <= y["high"] and y["low"] <= x["high"]:
   return True
  else:
   return False
 def intervalSearch(self,z):
  matches = self.intervalSearchAll(z)
  maxOverlapInt = self.nil
  maxOverlap = 0.0
  tieBreaker = 0.0
  y1 = z[0]
  y2 = z[1]
  for match in matches:
   x1 = match.key["low"]
   x2 = match.key["high"]
   if x1 == y1 and x2 == y2:
    return match
   overlap1 = (0.0 + min(x2,y2) - max(x1,y1)) / (x2-x1)
   overlap2 = (0.0 + min(x2,y2) - max(x1,y1)) / (y2-y1)
   overlap = max(overlap1,overlap2)
   minor = min(overlap1,overlap2)
   #if y2 == 4157430:
   # print str(x1) + " " + str(x2) + " " + str(overlap1) + " " + str(overlap2) + " " + str(overlap) + " " + str(maxOverlap) + " " + str(tieBreaker)
   if overlap == maxOverlap:
    if minor > tieBreaker:
     maxOverlap = overlap
     maxOverlapInt = match
     tieBreaker = minor
   elif overlap > maxOverlap:
    maxOverlap = overlap
    maxOverlapInt = match
    tieBreaker = minor
  return maxOverlapInt

 def intervalSearchAllRec(self,x,y,matches):
  if x == self.nil:
   return
  if self.overlap(y,x.key) == True:
   matches.append(x)
  if y["high"] < x.key["low"]:
   #print "Going Left"
   self.intervalSearchAllRec(x.left,y,matches)
  elif x.left != self.nil and y["low"] > x.left.max:
   #print "Going Right"
   self.intervalSearchAllRec(x.right,y,matches)
  else:
   #print "Going Both Ways"
   self.intervalSearchAllRec(x.left,y,matches)
   self.intervalSearchAllRec(x.right,y,matches)
   

 def intervalSearchAll(self,z):
  y = {"low":z[0],"high":z[1]}
  x = self.root
  matches = []
  self.intervalSearchAllRec(x,y,matches) 
  return matches

